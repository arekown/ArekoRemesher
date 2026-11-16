#!/usr/bin/env bash
# Builds the core and packages the Blender extension.
#   werkzeuge/bauen.sh                 build in org.freedesktop.Sdk//25.08 (matches Blender's runtime), run the check program,
#                                      copy the .so into the add-on, build the extension .zip into dist/
#   werkzeuge/bauen.sh --clang-libcxx  portability check: Clang + libc++ on the host (no packaging); the
#                                      toolchain root comes from ARM_CLANG_ROOT or from werkzeuge/lokal.sh
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SDK="org.freedesktop.Sdk//25.08"
ADDON="$REPO/addon/areko_remesher"

manifest_version() {
    sed -n 's/^version = "\(.*\)"$/\1/p' "$ADDON/blender_manifest.toml"
}

cmake_version() {
    sed -n 's/^project(ArekoRemesher VERSION \([0-9.]*\) .*/\1/p' "$REPO/CMakeLists.txt"
}

check_versions() {
    local m c
    m="$(manifest_version)"
    c="$(cmake_version)"
    if [[ -z "$m" || "$m" != "$c" ]]; then
        echo "version mismatch: blender_manifest.toml '$m', CMakeLists.txt '$c'" >&2
        exit 1
    fi
}

build_sdk() {
    flatpak run --filesystem="$REPO" --cwd="$REPO" --command=sh "$SDK" -c \
        'cmake --preset sdk && cmake --build --preset sdk && build/sdk/kern/arekoremesher_pruefung'
    mkdir -p "$ADDON/lib"
    # Copy beside it, then rename: a Blender that has the old library loaded keeps
    # its file and keeps running (overwriting the loaded file in place crashes it).
    cp "$REPO/build/sdk/kern/libarekoremesher.so" "$ADDON/lib/.libarekoremesher.so.neu"
    mv -f "$ADDON/lib/.libarekoremesher.so.neu" "$ADDON/lib/libarekoremesher.so"
}

package_extension() {
    mkdir -p "$REPO/dist"
    flatpak run --command=blender org.blender.Blender --factory-startup \
        --command extension build --source-dir "$ADDON" --output-dir "$REPO/dist"
}

build_clang_libcxx() {
    if [[ -z "${ARM_CLANG_ROOT:-}" && -f "$REPO/werkzeuge/lokal.sh" ]]; then
        source "$REPO/werkzeuge/lokal.sh"
    fi
    if [[ -z "${ARM_CLANG_ROOT:-}" ]]; then
        echo "set ARM_CLANG_ROOT to the root of a Clang toolchain with libc++" >&2
        exit 2
    fi
    export ARM_CLANG_ROOT
    (cd "$REPO" && cmake --preset clang-libcxx && cmake --build --preset clang-libcxx \
        && build/clang-libcxx/kern/arekoremesher_pruefung)
}

case "${1:-}" in
    "")
        check_versions
        build_sdk
        package_extension
        ;;
    --clang-libcxx)
        build_clang_libcxx
        ;;
    *)
        echo "usage: $0 [--clang-libcxx]" >&2
        exit 2
        ;;
esac
