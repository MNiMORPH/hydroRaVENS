#! /bin/sh

twine upload --repository testpypi dist/*
firefox https://test.pypi.org/project/hydroRaVENS/
