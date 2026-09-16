#!/bin/sh
# Optional local GAP 4.14.0 + CTblLib 1.3.9 installation.
set -eu
base_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
tool_dir="$base_dir/.tools"
mkdir -p "$tool_dir"
fetch() { [ -f "$tool_dir/$2" ] || curl -fL "$1" -o "$tool_dir/$2"; }
fetch https://github.com/gap-system/gap/releases/download/v4.14.0/gap-4.14.0-core.tar.gz gap-core.tar.gz
if [ ! -d "$tool_dir/gap-4.14.0" ]; then tar -xzf "$tool_dir/gap-core.tar.gz" -C "$tool_dir"; fi
mkdir -p "$tool_dir/gap-4.14.0/pkg"
fetch https://github.com/gap-system/gap/releases/download/v4.14.0/packages-required-v4.14.0.tar.gz gap-required.tar.gz
fetch https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/ctbllib-1.3.9.tar.gz ctbllib.tar.gz
fetch https://www.math.rwth-aachen.de/~Thomas.Breuer/atlasrep/atlasrep-2.1.9.tar.gz atlasrep.tar.gz
fetch https://github.com/gap-packages/utils/releases/download/v0.99/utils-0.99.tar.gz utils.tar.gz
for archive in gap-required.tar.gz ctbllib.tar.gz atlasrep.tar.gz utils.tar.gz; do
    tar -xzf "$tool_dir/$archive" -C "$tool_dir/gap-4.14.0/pkg"
done
if [ ! -x "$tool_dir/gap-4.14.0/gap" ]; then
    cd "$tool_dir/gap-4.14.0"
    if command -v brew >/dev/null 2>&1 && brew --prefix gmp >/dev/null 2>&1; then
        ./configure --with-gmp="$(brew --prefix gmp)"
    else
        ./configure
    fi
    make -j4
fi
printf '%s\n' "GAP ready at $tool_dir/gap-4.14.0/gap"
