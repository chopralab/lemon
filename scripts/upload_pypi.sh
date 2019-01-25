#!/bin/bash
echo -e "[distutils]" >> ~/.pypirc
echo -e "index-servers =" >> ~/.pypirc
echo -e "    test" >> ~/.pypirc
echo -e "" >> ~/.pypirc
echo -e "[test]" >> ~/.pypirc
echo -e "repository = https://test.pypi.org/legacy/" >> ~/.pypirc
echo -e "username = frodofine" >> ~/.pypirc
echo -e "password = $PYPI_PASSWORD" >> ~/.pypirc

pip install twine
twine upload --repository test dist/*
