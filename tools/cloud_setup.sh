#!/bin/bash
## Environment setup for PSMFC-mytilus-byssus-pilot on Ubuntu 24.04 ("noble"), written as the
## setup script of the Claude Code cloud environment (paste it into the environment's settings,
## Setup script; new sessions run it) and usable on any Ubuntu 24.04 machine with root:
##   1. R 4.6.1 from CRAN's Ubuntu repository, with R's recommended packages built for R 4.6;
##   2. the system libraries the R packages compile against, and pandoc;
##   3. BLAST+ 2.15.0 from NCBI (checksum-verified), on the PATH;
##   4. the R packages of renv.lock (R 4.6.1, Bioconductor 3.23) in their own library,
##      /opt/R/site-library-4.6, which R uses instead of the image's R 4.3 packages.
## Safe to run again: each part is skipped or brought up to date if already done.
## Step 4 builds about 280 packages from source, about an hour on 4 cores; set
## RESTORE_R_PACKAGES=0 to skip it and run it later in the session (the command is printed).
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive

R_VERSION=4.6.1
R_LIB=/opt/R/site-library-4.6
BLAST_VERSION=2.15.0
BLAST_MD5=0abb189643afd79f3fbbd2bf27db6428        # ncbi-blast-2.15.0+-x64-linux.tar.gz (NCBI .md5)
## the repository: REPO_DIR if set, else the checkout this script sits in (tools/..), else the
## cloud environment's checkout
here=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" 2>/dev/null && pwd || true)
if [ -z "${REPO_DIR:-}" ]; then
  if [ -n "${here}" ] && [ -f "${here}/../renv.lock" ]; then REPO_DIR=$(cd "${here}/.." && pwd)
  else REPO_DIR=/home/user/PSMFC-mytilus-byssus-pilot; fi
fi
LOCK_BRANCH=${LOCK_BRANCH:-claude/go-db-2026}     # renv.lock is on this branch until PR 3 is merged
RESTORE_R_PACKAGES=${RESTORE_R_PACKAGES:-1}
SUDO=""; [ "$(id -u)" -eq 0 ] || SUDO=sudo
say() { printf '\n== %s\n' "$*"; }

say "1. R ${R_VERSION} from CRAN's Ubuntu repository"
$SUDO apt-get update -qq
$SUDO apt-get install -y -qq --no-install-recommends ca-certificates curl gnupg git
curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
  | $SUDO gpg --dearmor --yes -o /usr/share/keyrings/cran-r.gpg
echo "deb [signed-by=/usr/share/keyrings/cran-r.gpg] https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" \
  | $SUDO tee /etc/apt/sources.list.d/cran-r.list > /dev/null
$SUDO apt-get update -qq
if ! $SUDO apt-get install -y -qq --no-install-recommends "r-base-core=${R_VERSION}-*" "r-base-dev=${R_VERSION}-*"; then
  echo "R ${R_VERSION} is no longer offered; installing the newest R 4.6 instead" >&2
  $SUDO apt-get install -y -qq --no-install-recommends r-base-core r-base-dev
fi
## R's recommended packages (Matrix, MASS, mgcv, survival, ...): the image has R 4.3 builds
$SUDO apt-get install -y -qq --no-install-recommends r-recommended r-cran-boot r-cran-class \
  r-cran-cluster r-cran-codetools r-cran-foreign r-cran-kernsmooth r-cran-lattice r-cran-mass \
  r-cran-matrix r-cran-mgcv r-cran-nlme r-cran-nnet r-cran-rpart r-cran-spatial r-cran-survival

say "2. System libraries for the R packages"
$SUDO apt-get install -y -qq --no-install-recommends \
  libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype-dev libpng-dev libtiff-dev libjpeg-dev libwebp-dev libcairo2-dev libxt-dev \
  libglpk-dev libgmp-dev libicu-dev libmagick++-dev libuv1-dev libnlopt-dev libgit2-dev \
  libbz2-dev liblzma-dev zlib1g-dev libdeflate-dev pandoc

say "3. BLAST+ ${BLAST_VERSION}"
BLAST_DIR=/opt/ncbi-blast-${BLAST_VERSION}+
if ! "${BLAST_DIR}/bin/blastx" -version 2>/dev/null | grep -q "blastx: ${BLAST_VERSION}+"; then
  tmp=$(mktemp -d); tarball=ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
  curl -fsS -o "${tmp}/${tarball}" "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/${tarball}"
  echo "${BLAST_MD5}  ${tmp}/${tarball}" | md5sum -c --quiet
  $SUDO tar -xzf "${tmp}/${tarball}" -C /opt
  rm -rf "${tmp}"
fi
for b in "${BLAST_DIR}"/bin/*; do $SUDO ln -sf "$b" /usr/local/bin/; done
blastx -version | head -1

say "4. R packages of renv.lock in ${R_LIB}"
$SUDO mkdir -p "${R_LIB}"; $SUDO chmod 0775 "${R_LIB}"
## R looks only in this library and its own (not /usr/lib/R/site-library, whose R 4.3 builds
## fail to load in R 4.6). This line overrides an R_LIBS_SITE set in the environment; put
## another library first with R_LIBS.
grep -qx "R_LIBS_SITE=\"${R_LIB}\"" /etc/R/Renviron.site \
  || echo "R_LIBS_SITE=\"${R_LIB}\"" | $SUDO tee -a /etc/R/Renviron.site > /dev/null
lock=""
if [ -f "${REPO_DIR}/renv.lock" ]; then
  lock="${REPO_DIR}/renv.lock"
elif git -C "${REPO_DIR}" fetch -q origin "${LOCK_BRANCH}" 2>/dev/null; then
  lock=$(mktemp --suffix=.lock); git -C "${REPO_DIR}" show FETCH_HEAD:renv.lock > "${lock}"
fi
restore_cmd="MAKEFLAGS=-j$(nproc) Rscript -e 'install.packages(\"renv\", repos = \"https://cloud.r-project.org\", lib = \"${R_LIB}\"); renv::restore(lockfile = \"<renv.lock>\", library = \"${R_LIB}\", prompt = FALSE)'"
if [ -z "${lock}" ]; then
  echo "No renv.lock found (looked in ${REPO_DIR} and on ${LOCK_BRANCH}). Later, in the repository:" >&2
  echo "  ${restore_cmd}" >&2
elif [ "${RESTORE_R_PACKAGES}" != "1" ]; then
  echo "Skipped (RESTORE_R_PACKAGES=${RESTORE_R_PACKAGES}). Later, in the repository:"
  echo "  ${restore_cmd}"
else
  MAKEFLAGS="-j$(nproc)" $SUDO Rscript -e "
    install.packages('renv', repos = 'https://cloud.r-project.org', lib = '${R_LIB}', quiet = TRUE)
    renv::restore(lockfile = '${lock}', library = '${R_LIB}', prompt = FALSE)"
  Rscript -e "
    lock <- jsonlite::fromJSON('${lock}', simplifyVector = FALSE)
    want <- vapply(lock\$Packages, \`[[\`, '', 'Version')
    have <- installed.packages()[, 'Version']
    bad <- names(want)[is.na(have[names(want)]) | have[names(want)] != want]
    cat(R.version.string, '| Bioconductor', as.character(BiocManager::version()), '|', length(want), 'packages,',
        if (length(bad)) paste('NOT as in renv.lock:', paste(bad, collapse = ', ')) else 'all as in renv.lock', '\n')
    if (length(bad)) quit(status = 1)"
fi
say "Done"
