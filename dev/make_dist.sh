#!/bin/bash

set -e

GAP_PKG_RELEASE_DATE=$(date -I) ./dev/release-gap-package --skip-existing-release --release-script dev/.release
