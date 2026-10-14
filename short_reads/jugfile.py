"""
Runs fargene/rgi/deeparg/resfinder over the samples listed in data/samples.txt.

    pixi run jug-status
    pixi run jug-execute

Safe to re-run: jug skips tasks it already finished, even if you add new
samples in between runs. Only settings that are passed to the tasks as
arguments (the directories, CARD_JSON, CARD_VERSION, RESFINDER_DB_VERSION,
RESFINDER_MIN_COV, RESFINDER_THRESHOLD) are part of the task hashes; changing
THREADS, TMPDIR or FARGENE_MODELS does *not* invalidate anything, so use `jug
invalidate` if the results should be redone.

Each tool lives in its own pixi environment (see pixi.toml) with its own
binary -- this process (the 'jug' environment) only needs jug itself and
shells out via `pixi run -e <env> ...` for the actual work.

The one-time database setup is part of the graph too: prepare_rgi_db builds
CARD's localDB, prepare_deeparg_db downloads the DeepARG model bundle, and
prepare_resfinder_db clones and indexes resfinder_db, each once, with the
run_* tasks that need them waiting on the result.

Every tool invocation gets a scratch directory of its own (under
TMPDIR), laid out as

    <scratch>/reads/preproc.pair.{1,2}.fq[.gz]   quality-trimmed reads
    <scratch>/out/                               what the tool is told to write
    <scratch>/localDB -> rgi_db/localDB          (rgi only)

The reads are produced there by running preprocess.ngl through ngless, the
tool then runs with the scratch directory as its working directory, and
finally <scratch>/out is copied to OUTPUT_DIR/<sample>/<tool>/ and the whole
scratch directory is deleted. Everything a tool leaves lying around outside
that out/ directory -- deeparg's trimmomatic/vsearch intermediates, rgi's
bowtie2 scratch -- therefore dies with the scratch directory instead of
accumulating in the output tree, and OUTPUT_DIR/<sample>/<tool>/ holds
nothing but that tool's results.

Only the files listed in KEEP are published, gzipped (see publish_results);
everything else the tools write, inside out/ or not, dies with the scratch
directory. fargene is additionally pointed at a tmp dir in scratch and rgi
runs with --clean, so that their intermediates don't pile up in scratch while
the tool is still running.
"""
import contextlib
import fnmatch
import gzip
import os
import shutil
import re
import subprocess
import tempfile

from jug import TaskGenerator

# Old-fashioned NGLess-style organisation: data/samples.txt lists one sample
# name per line, and each sample's reads are in METAGENOMES_DIR/<name>/, where
# ngless' load_fastq_directory picks them up (.fq/.fastq, optionally
# compressed, paired up on a .1/.2 or _1/_2 suffix), e.g.
#     data/metagenomes/s1/s1.pair.1.fq.gz
#     data/metagenomes/s1/s1.pair.2.fq.gz
SAMPLES = [line.strip() for line in open("data/samples.txt")]
METAGENOMES_DIR = "data/metagenomes"

OUTPUT_DIR = "output"

# Parent directory for the per-invocation scratch directories, which hold the
# preprocessed reads plus whatever the tool scribbles next to them. None means
# tempfile's default ($TMPDIR, else /tmp). On a cluster, point this at
# node-local scratch: every byte written here is thrown away once the tool's
# results have been copied into OUTPUT_DIR.
TMPDIR = None

THREADS = 4

# `rgi --local` resolves its database as localDB/ relative to the directory it
# is run from; prepare_rgi_db builds it once under RGI_LOCALDB_DIR and run_rgi
# symlinks it into each run's scratch directory. Gitignored.
RGI_LOCALDB_DIR = "rgi_db"
# Path to a local card.json to build the database from; None downloads CARD
# release CARD_VERSION from card.mcmaster.ca instead.
CARD_JSON = None
CARD_VERSION = "4.0.0"
# Where prepare_deeparg_db puts the DeepARG v2 model+database bundle. Gitignored.
DEEPARG_HF_DIR = "deeparg_hf"
# Where prepare_resfinder_db clones and KMA-indexes resfinder_db. Gitignored.
RESFINDER_DB_DIR = "resfinder_db"
RESFINDER_DB_URL = "https://bitbucket.org/genomicepidemiology/resfinder_db.git"
# A tag of resfinder_db; 2.4.0 is the version used for the GMGC10 runs (and
# the one bundled with ResFinder 4.6.0, the version in pixi.toml).
RESFINDER_DB_VERSION = "2.4.0"
RESFINDER_DB_STAMP = ".install-complete"
# Minimum breadth of coverage (-l) and identity (-t), both as fractions: the
# same values as in resfinder_dna/Snakefile_resfinder_dna.
RESFINDER_MIN_COV = 0.6
RESFINDER_THRESHOLD = 0.8

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PIXI_MANIFEST = os.path.join(SCRIPT_DIR, "pixi.toml")
PREPROCESS_NGL = os.path.join(SCRIPT_DIR, "preprocess.ngl")
DOWNLOAD_DEEPARG_DB = os.path.join(SCRIPT_DIR, "download_deeparg_db.py")


# Only for the stamp-file name; the script's huggingface_hub import lives
# inside its own function, so importing it here (under 'jug') is fine. It is
# *run* below via `pixi run -e deeparg`, never called from this process.
from download_deeparg_db import STAMP_FILE as DEEPARG_DB_STAMP

FARGENE_MODELS = {
    "class_a": "class_a", "class_b1_b2": "class_b_1_2", "class_b3": "class_b_3",
    "class_c": "class_c", "class_d1": "class_d_1", "class_d2": "class_d_2",
    "mph": "mph", "erm_1": "erm_type_a", "erm_2": "erm_type_f",
    "tet_enzyme": "tet_enzyme", "tet_rpg": "tet_rpg",
    "aph2b": "aminoglycoside_model_g", "aph3p": "aminoglycoside_model_h", "aph6": "aminoglycoside_model_i",
    "aac2p": "aminoglycoside_model_a", "aac3_1": "aminoglycoside_model_b", "aac3_2": "aminoglycoside_model_c",
    "aac6p_1": "aminoglycoside_model_d", "aac6p_2": "aminoglycoside_model_e", "aac6p_3": "aminoglycoside_model_f",
    "tet_efflux": "tet_efflux", "qnr": "qnr",
}

# Basename handed to preprocess.ngl; ngless turns it into <base>.pair.1.<ext>,
# <base>.pair.2.<ext> and (only for samples that carry unpaired reads)
# <base>.singles.<ext>.
PREPROC_BASE = "preproc"


def pixi_run(env, args, **kwargs):
    """Run `args` inside the `env` pixi environment, installing it if needed.

    This process runs under the 'jug' environment, so it can't see the tools;
    `pixi run -e <env>` activates theirs, which also lets fargene/rgi/deeparg
    find the sibling binaries they shell out to (deeparg -> trimmomatic, rgi
    -> bamtools). --manifest-path is what makes `cwd` free to be anything:
    without it pixi locates pixi.toml by searching upwards from the working
    directory, which would rule out the scratch directories below (they live
    outside the workspace, typically on node-local scratch). The commands are
    passed as plain commands rather than as pixi tasks for the same reason:
    tasks always run from the workspace root.
    """
    kwargs.setdefault("check", True)
    return subprocess.run(["pixi", "run", "--manifest-path", PIXI_MANIFEST, "-e", env, *args],
                          **kwargs)


# The files of each tool's out/ directory that are published, as fnmatch
# patterns on the path relative to out/ (where `*` also matches `/`). They
# name the files as the tools write them; they are gzipped on publishing, so
# e.g. sample.bwt.gene_mapping_data.txt ends up as .txt.gz.
KEEP = {
    "rgi": [
        "sample.bwt.gene_mapping_data.txt",
        "sample.bwt.allele_mapping_data.txt",
    ],
    "deeparg": [
        "sample.clean.deeparg.align.daa.tsv",
        "sample.clean.deeparg.mapping.ARG",
        "sample.clean.deeparg.mapping.ARG.merged",
        "sample.clean.deeparg.mapping.ARG.merged.quant",
    ],
    # Per model (<class>/): the reconstructed genes, and the reads passing
    # the HMM -- per mate, after trim_galore, and concatenated over inputs.
    "fargene": [
        "*/predictedGenes/predicted-orfs.fasta",
        "*/retrievedFragments/*_1_retrieved.fastq",
        "*/retrievedFragments/*_2_retrieved.fastq",
        "*/retrievedFragments/trimmedReads/*_1_retrieved_val_1.fq",
        "*/retrievedFragments/trimmedReads/*_2_retrieved_val_2.fq",
        "*/retrievedFragments/all_retrieved_1.fastq",
        "*/retrievedFragments/all_retrieved_2.fastq",
    ],
    # The hits (flat, and split up by drug class), the sequence of the reads'
    # consensus over each hit, the JSON with everything (per-hit depth
    # included), and KMA's unfiltered per-template results, which carry the
    # depth too.
    "resfinder": [
        "ResFinder_results_tab.txt",
        "ResFinder_results_table.txt",
        "ResFinder_Hit_in_genome_seq.fsa",
        "sample.json",
        "resfinder_kma/kma_*.res",
    ],
}


def trim_results(out_dir, keep):
    """Delete every file under `out_dir` not matching one of the `keep`
    patterns (see KEEP), then any directories left empty."""
    for root, dirs, files in os.walk(out_dir, topdown=False):
        for fn in files:
            path = os.path.join(root, fn)
            # A .gz is matched by its uncompressed name, so that this also
            # applies to already-published results.
            rel = os.path.relpath(path, out_dir).removesuffix(".gz")
            if not any(fnmatch.fnmatch(rel, pat) for pat in keep):
                os.remove(path)
        for d in dirs:
            path = os.path.join(root, d)
            if not os.path.islink(path) and not os.listdir(path):
                os.rmdir(path)


def compress_results(out_dir):
    """Gzip every file under `out_dir` in place (foo -> foo.gz)."""
    for root, _, files in os.walk(out_dir):
        for fn in files:
            path = os.path.join(root, fn)
            if fn.endswith(".gz"):
                continue
            with open(path, "rb") as src, \
                    gzip.open(path + ".gz", "wb", compresslevel=6) as dst:
                shutil.copyfileobj(src, dst, 16 * 1024 * 1024)
            os.remove(path)


def publish_results(src, dest, keep):
    """Trim and compress the scratch directory's out/ tree and copy it to
    OUTPUT_DIR.

    Staged next to `dest` and renamed into place so that an interrupted copy
    never leaves a half-written result looking like a finished one, and so
    that re-running an invalidated task replaces the old results wholesale
    rather than merging into them.
    """
    trim_results(src, keep)
    compress_results(src)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    staging = f"{dest}.incoming.{os.getpid()}"
    superseded = f"{dest}.superseded.{os.getpid()}"
    shutil.rmtree(staging, ignore_errors=True)
    shutil.copytree(src, staging)
    if os.path.exists(dest):
        os.rename(dest, superseded)
    os.rename(staging, dest)
    shutil.rmtree(superseded, ignore_errors=True)
    return dest


@contextlib.contextmanager
def tool_scratch(sample_dir, tool, compressed):
    """Preprocess one sample into a fresh scratch directory, for one tool run.

    Yields (scratch, out_dir, r1, r2): the directory the tool should run in,
    the directory it should be told to write to (whose KEEP files are copied
    out by publish_results), and the two mate files. `compressed` picks the
    output extension -- fargene cannot read gzipped FastQ, so it asks for
    plain .fq and lets ngless decompress on the way out.

    The whole tree is removed on the way out, including on failure: jug runs
    without --keep-failed, so a failed task is retried from scratch anyway and
    there is nothing here worth keeping. Cleanup errors are deliberately not
    suppressed: failing to remove the tree usually means something else went
    wrong (a tool still holding files open, a read-only file it left behind).
    """
    if not os.path.isdir(sample_dir):
        raise RuntimeError(f"no read directory for this sample: {sample_dir}")
    # fargene quality-filters its input with trim_galore, and when that fails
    # (as it does on some very old systems) fargene carries on regardless and
    # produces empty results. Fail loudly up front instead.
    pixi_run('fargene', [
        'trim_galore', '--version'])
    tmp_parent = os.path.abspath(TMPDIR) if TMPDIR else None
    if tmp_parent:
        os.makedirs(tmp_parent, exist_ok=True)
    with tempfile.TemporaryDirectory(
            prefix=f"{tool}.{os.path.basename(sample_dir)}.",
            dir=tmp_parent) as scratch:
        reads_dir = os.path.join(scratch, "reads")
        out_dir = os.path.join(scratch, "out")
        os.makedirs(reads_dir)
        os.makedirs(out_dir)

        ext = ".fq.gz" if compressed else ".fq"
        pixi_run("ngless", [
            "ngless",
            "--jobs", str(THREADS),
            # Keep ngless' own scratch and its QC report inside ours (the
            # default would be a single preprocess.ngl.output_ngless/ next to
            # the script, which concurrent invocations would fight over).
            "--temporary-directory", scratch,
            "-o", os.path.join(scratch, "ngless-report"),
            PREPROCESS_NGL,
            sample_dir,
            os.path.join(reads_dir, PREPROC_BASE + ext),
        ], cwd=scratch)

        r1 = os.path.join(reads_dir, f"{PREPROC_BASE}.pair.1{ext}")
        r2 = os.path.join(reads_dir, f"{PREPROC_BASE}.pair.2{ext}")
        missing = [p for p in (r1, r2) if not os.path.isfile(p)]
        if missing:
            # ngless omits the paired outputs entirely for a single-end
            # sample, which none of these three tools are wired up for here.
            raise RuntimeError(
                f"preprocessing {sample_dir} did not produce paired reads "
                f"(missing {', '.join(os.path.basename(p) for p in missing)}); "
                f"ngless wrote {sorted(os.listdir(reads_dir))}")

        yield scratch, out_dir, r1, r2


@TaskGenerator
def run_fargene(sample_dir, outdir):
    with tool_scratch(sample_dir, "fargene", compressed=False) as (scratch, out_dir, r1, r2):
        for class_name, model in FARGENE_MODELS.items():
            # By default fargene keeps its intermediates (including a FASTA
            # copy of both read files, >500M per model per sample) in
            # <-o>/tmpdir; keep them in scratch and drop them per model.
            model_tmp = os.path.join(scratch, "fargene-tmp", class_name)
            pixi_run("fargene", [
                "fargene",
                "-i", r1, r2,
                "--hmm-model", model,
                "--meta",
                "-o", os.path.join(out_dir, class_name),
                "--tmp-dir", model_tmp,
                "-p", str(THREADS),
                "--force",
            ], cwd=scratch)
            shutil.rmtree(model_tmp)
        return publish_results(out_dir, os.path.join(outdir, "fargene"),
                               KEEP["fargene"])


@TaskGenerator
def prepare_rgi_db(work_dir, card_json, card_version):
    # Returns work_dir itself (not work_dir/localDB) -- run_rgi symlinks
    # work_dir/localDB into its own scratch directory, and the build steps
    # below likewise run one level above localDB/, which is where `rgi
    # ... --local` looks for it.
    os.makedirs(work_dir, exist_ok=True)
    localdb = os.path.join(work_dir, "localDB")
    if os.path.isdir(localdb):
        # Already built -- nothing to do.
        return work_dir

    if card_json:
        if not os.path.isfile(card_json):
            raise RuntimeError(f"card.json not found: {card_json}")
        shutil.copyfile(card_json, os.path.join(work_dir, "card.json"))
    else:
        url = f"https://card.mcmaster.ca/download/0/broadstreet-v{card_version}.tar.bz2"
        tarball = os.path.join(work_dir, "card_data.tar.bz2")
        subprocess.run(["curl", "-fsSL", "-o", tarball, url], check=True)
        subprocess.run(["tar", "xjf", os.path.basename(tarball), "./card.json"], cwd=work_dir, check=True)
        os.remove(tarball)

    log_path = os.path.join(work_dir, "card_annotation.log")
    with open(log_path, "wb") as log:
        pixi_run("rgi", ["rgi", "card_annotation", "-i", "card.json"],
                 cwd=work_dir, stdout=log, stderr=subprocess.STDOUT)

    pattern = re.compile(r"^card_database_v(.+)\.fasta$")
    versions = sorted({
        m.group(1) for fn in os.listdir(work_dir)
        if (m := pattern.match(fn)) and "_all" not in fn
    })
    if len(versions) != 1:
        raise RuntimeError(
            f"expected exactly one card_database_v*.fasta (excluding _all), "
            f"found {len(versions)} -- see {log_path}")
    ver = versions[0]

    pixi_run("rgi", ["rgi", "clean", "--local"], cwd=work_dir)
    pixi_run("rgi", [
        "rgi", "load",
        "--card_json", "card.json",
        "--card_annotation", f"card_database_v{ver}.fasta",
        "--card_annotation_all_models", f"card_database_v{ver}_all.fasta",
        "--local",
    ], cwd=work_dir)

    seed_fna = os.path.join(work_dir, "_seed.fna")
    with open(seed_fna, "w") as fh:
        fh.write(">seed\nATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT\n")
    pixi_run("rgi", [
        "rgi", "main", "-a", "DIAMOND", "-i", "_seed.fna",
        "-o", "_seed_out", "--local", "--clean", "-t", "contig", "-n", "1",
    ], cwd=work_dir)
    for fn in os.listdir(work_dir):
        if fn.startswith("_seed"):
            p = os.path.join(work_dir, fn)
            shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)

    return work_dir


@TaskGenerator
def run_rgi(sample_dir, outdir, rgi_work_dir):
    with tool_scratch(sample_dir, "rgi", compressed=True) as (scratch, out_dir, r1, r2):
        # `rgi ... --local` resolves its database as localDB/ under the
        # working directory, so point one there. A symlink (rather than a
        # copy) keeps this cheap and lets concurrent samples share the build.
        os.symlink(os.path.join(rgi_work_dir, "localDB"),
                   os.path.join(scratch, "localDB"))
        pixi_run("rgi", [
            "rgi", "bwt",
            "--read_one", r1,
            "--read_two", r2,
            "--output_file", os.path.join(out_dir, "sample.bwt"),
            "--local",
            "--threads", str(THREADS),
            # Removes every <output_file>*.temp* file: the SAM, the unsorted
            # and sorted-unfiltered BAMs, the per-read tables, ... (>1G per
            # sample) as soon as rgi is done with them.
            "--clean",
        ], cwd=scratch)
        return publish_results(out_dir, os.path.join(outdir, "rgi"),
                               KEEP["rgi"])


@TaskGenerator
def prepare_deeparg_db(hf_dir):
    """Fetch the DeepARG model+database bundle into `hf_dir`; returns it.

    The download is a separate script run under the 'deeparg' environment
    because it needs huggingface_hub, which this process (the 'jug'
    environment) doesn't have -- see download_deeparg_db.py, which `pixi run
    download-deeparg-db <dir>` runs by hand in exactly the same way.

    A directory carrying the stamp file is taken as already downloaded and
    left alone, so a bundle staged by hand (or by an earlier, since
    invalidated, run of this task) is not re-fetched.
    """
    if os.path.isfile(os.path.join(hf_dir, DEEPARG_DB_STAMP)):
        return hf_dir

    pixi_run("deeparg", ["python", DOWNLOAD_DEEPARG_DB, hf_dir])
    if not os.path.isfile(os.path.join(hf_dir, DEEPARG_DB_STAMP)):
        raise RuntimeError(
            f"{DOWNLOAD_DEEPARG_DB} exited successfully but left no "
            f"{DEEPARG_DB_STAMP} in {hf_dir}")
    return hf_dir


@TaskGenerator
def run_deeparg(sample_dir, outdir, hf_dir):
    with tool_scratch(sample_dir, "deeparg", compressed=True) as (scratch, out_dir, r1, r2):
        pixi_run("deeparg", [
            "deeparg", "short_reads_pipeline",
            "--forward_pe_file", r1,
            "--reverse_pe_file", r2,
            "--output_file", os.path.join(out_dir, "sample"),
            # Without this deeparg downloads the bundle itself, once per run;
            # prepare_deeparg_db has already put it here.
            "--hf-model-path", hf_dir,
        ], cwd=scratch)
        return publish_results(out_dir, os.path.join(outdir, "deeparg"),
                               KEEP["deeparg"])


@TaskGenerator
def prepare_resfinder_db(db_dir, db_version):
    """Clone resfinder_db at tag `db_version` into `db_dir` and index it with
    KMA (which is what ResFinder uses on reads); returns `db_dir`.

    Only resfinder_db is needed: run_resfinder asks for acquired genes only,
    so neither pointfinder_db nor disinfinder_db is ever read.

    A directory carrying the stamp file is taken as already built and left
    alone; one without it is a leftover of a failed build and is redone.
    """
    stamp = os.path.join(db_dir, RESFINDER_DB_STAMP)
    if os.path.isfile(stamp):
        return db_dir
    shutil.rmtree(db_dir, ignore_errors=True)
    pixi_run("resfinder", [
        "git", "-c", "advice.detachedHead=false",
        "clone", "--quiet", "--depth", "1", "--branch", db_version,
        RESFINDER_DB_URL, db_dir])
    with open(os.path.join(db_dir, "VERSION")) as fh:
        found = fh.read().strip()
    if found != db_version:
        raise RuntimeError(
            f"resfinder_db tag {db_version} says it is version {found}")

    # This is what the database's own INSTALL.py does (one index per drug
    # class listed in config, plus all.fsa). `kma index` (1.6.17, at least)
    # exits with a stale errno -- 17, EEXIST -- even when it succeeds, which
    # is also why INSTALL.py's os.system ignores its status; check that the
    # index files are there instead.
    with open(os.path.join(db_dir, "config")) as fh:
        drugs = [line.split("\t")[0].strip() for line in fh
                 if line.strip() and not line.startswith("#")]
    for drug in drugs + ["all"]:
        pixi_run("resfinder", ["kma", "index", "-i", f"{drug}.fsa", "-o", drug],
                 cwd=db_dir, stdout=subprocess.DEVNULL, check=False)
        missing = [f"{drug}{ext}" for ext in (".comp.b", ".length.b", ".name", ".seq.b")
                   if not os.path.isfile(os.path.join(db_dir, f"{drug}{ext}"))]
        if missing:
            raise RuntimeError(f"`kma index` did not produce {', '.join(missing)}")

    open(stamp, "w").close()
    return db_dir


@TaskGenerator
def run_resfinder(sample_dir, outdir, db_dir, min_cov, threshold):
    with tool_scratch(sample_dir, "resfinder", compressed=True) as (scratch, out_dir, r1, r2):
        pixi_run("resfinder", [
            "python", "-m", "resfinder",
            "--inputfastq", r1, r2,
            "--outputPath", out_dir,
            # Otherwise named after the first read file, i.e. preproc.json.
            "--out_json", os.path.join(out_dir, "sample.json"),
            "--acquired",
            "--db_path_res", db_dir,
            "--min_cov", str(min_cov),
            "--threshold", str(threshold),
        ], cwd=scratch)
        return publish_results(out_dir, os.path.join(outdir, "resfinder"),
                               KEEP["resfinder"])


# Prepared once, ahead of the samples; every run_rgi/run_deeparg/run_resfinder
# task takes the corresponding directory as an argument and so waits for it.
rgi_work_dir = prepare_rgi_db(os.path.abspath(RGI_LOCALDB_DIR),
                              CARD_JSON, CARD_VERSION)
deeparg_hf_dir = prepare_deeparg_db(os.path.abspath(DEEPARG_HF_DIR))
resfinder_db_dir = prepare_resfinder_db(os.path.abspath(RESFINDER_DB_DIR),
                                        RESFINDER_DB_VERSION)

metagenomes_dir = os.path.abspath(METAGENOMES_DIR)
output_dir = os.path.abspath(OUTPUT_DIR)

# Adding another tool is a run_<tool> task above plus one line here: jug
# schedules the new tasks and leaves every result already on disk untouched.
for sample in SAMPLES:
    # Absolute, and passed as task arguments rather than looked up inside the
    # tasks, so that moving METAGENOMES_DIR invalidates the affected tasks.
    sample_dir = os.path.join(metagenomes_dir, sample)
    outdir = os.path.join(output_dir, sample)

    run_fargene(sample_dir, outdir)
    run_rgi(sample_dir, outdir, rgi_work_dir)
    run_deeparg(sample_dir, outdir, deeparg_hf_dir)
    run_resfinder(sample_dir, outdir, resfinder_db_dir,
                  RESFINDER_MIN_COV, RESFINDER_THRESHOLD)
