#!/usr/bin/env bash
# Builds the core and packages the Blender extension, one .zip per platform.
#   werkzeuge/bauen.sh                 build in org.freedesktop.Sdk//25.08 (matches Blender's runtime), run the check program,
#                                      copy the .so into the add-on, build dist/areko_remesher-<version>-linux-x64.zip
#   werkzeuge/bauen.sh --other-platforms [run-id]
#                                      download the Windows and macOS libraries built by GitHub Actions
#                                      (.github/workflows/kern.yml; default: the latest successful run for HEAD)
#                                      and build dist/areko_remesher-<version>-windows-x64.zip and -macos.zip
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

# package_platform <suffix> <platforms as TOML list> <library file>
# Copies the add-on into build/paket/<suffix>, keeps only the given library in
# lib/ and only the given platforms in the manifest, then builds the .zip.
package_platform() {
    local suffix="$1" platforms="$2" library="$3"
    local stage="$REPO/build/paket/$suffix/areko_remesher"
    rm -rf "$REPO/build/paket/$suffix"
    mkdir -p "$stage"
    cp -r "$ADDON/." "$stage/"
    rm -rf "$stage/lib"
    mkdir "$stage/lib"
    cp "$library" "$stage/lib/"
    sed -i "s/^platforms = .*/platforms = $platforms/" "$stage/blender_manifest.toml"
    grep -qxF "platforms = $platforms" "$stage/blender_manifest.toml"
    mkdir -p "$REPO/dist"
    flatpak run --command=blender org.blender.Blender --factory-startup \
        --command extension build --source-dir "$stage" \
        --output-filepath "$REPO/dist/areko_remesher-$(manifest_version)-$suffix.zip"
}

package_other_platforms() {
    local run="${1:-}"
    if [[ -z "$run" ]]; then
        run="$(cd "$REPO" && gh run list --workflow kern.yml --commit "$(git rev-parse HEAD)" \
            --status success --limit 1 --json databaseId --jq '.[0].databaseId')"
        if [[ -z "$run" ]]; then
            echo "no successful run of kern.yml for HEAD" >&2
            exit 1
        fi
    fi
    local artifacts="$REPO/build/artefakte"
    rm -rf "$artifacts"
    mkdir -p "$artifacts"
    (cd "$REPO" && gh run download "$run" --dir "$artifacts")
    package_platform windows-x64 '["windows-x64"]' "$artifacts/arekoremesher-windows-x64/arekoremesher.dll"
    package_platform macos '["macos-arm64", "macos-x64"]' "$artifacts/arekoremesher-macos-universal/libarekoremesher.dylib"
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
        package_platform linux-x64 '["linux-x64"]' "$ADDON/lib/libarekoremesher.so"
        ;;
    --other-platforms)
        check_versions
        package_other_platforms "${2:-}"
        ;;
    --clang-libcxx)
        build_clang_libcxx
        ;;
    *)
        echo "usage: $0 [--clang-libcxx | --other-platforms [run-id]]" >&2
        exit 2
        ;;
esac
