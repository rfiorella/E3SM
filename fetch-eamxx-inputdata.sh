#!/usr/bin/env bash
# Download the EAMxx input files listed in eamxx-input-files.txt into /inputdata,
# preserving their relative paths. Intended to be called from a Dockerfile RUN
# layer, e.g.:
#
#   COPY eamxx-input-files.txt fetch-eamxx-inputdata.sh /tmp/
#   RUN /tmp/fetch-eamxx-inputdata.sh
set -euo pipefail

BASE_URL="https://web.lcrc.anl.gov/public/e3sm/inputdata/atm"
FILE_LIST="$(dirname "$(readlink -f "$0")")/eamxx-input-files.txt"
DEST_ROOT="/inputdata"

mkdir -p "$DEST_ROOT"

while IFS= read -r rel_path || [ -n "$rel_path" ]; do
  [ -z "$rel_path" ] && continue

  dest_path="${DEST_ROOT}/${rel_path}"
  mkdir -p "$(dirname "$dest_path")"

  wget --no-verbose --timestamping --directory-prefix="$(dirname "$dest_path")" \
    "${BASE_URL}/${rel_path}"
done < "$FILE_LIST"
