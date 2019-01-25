#!/bin/bash
pip install twine
SCRIPT_DIR=$(cd $(dirname $0) || exit 1; pwd)
twine upload --repository-url https://test.pypi.org/ dist/*