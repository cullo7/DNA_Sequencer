HOST=127.0.0.1
TEST_PATH=./

clean-pyc:
	rm -rf '*.pyc'
	rm -rf '*.pyo'
	rm -rf '*~'

clean-build:
	rm -rf __pycache__/

sort:
	sh -c "isort --skip-glob=.tox --recursive . "

lint:
	flake8 --exclude=.tox

run:
	python3 compiler.py

irun:
	python3 compiler.py < input.txt

.PHONY: clean-pyc clean-build

help:
	@echo "    clean-pyc"
	@echo "        Remove python artifacts."
	@echo "    clean-build"
	@echo "        Remove build artifacts."
	@echo "    isort"
	@echo "        Sort import statements."
	@echo "    lint"
	@echo "        Check style with flake8."
	@echo '    run'
	@echo '        Run the `my_project` service on your local machine.'

