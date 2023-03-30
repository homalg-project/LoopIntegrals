#!/bin/bash

set -e

GAP_PKG_RELEASE_DATE=$(date -I) ./dev/release-gap-package --release-script dev/.release --only-tarball
