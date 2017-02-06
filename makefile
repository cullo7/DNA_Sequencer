HOST=127.0.0.1
TEST_PATH=./

run:
	python3 main.py

clean:
	make clean-pyc && make clean-build

clean-pyc:
	rm -rf *.pyc
	rm -rf *.pyo
	rm -rf *~
	rm -rf output/*

clean-build:
	rm -rf scripts/_pycache__/

sort:
	sh -c "isort --skip-glob=.tox --recursive . "

lint:
	flake8 --exclude=.tox

irun:
	python3 main.py < input.txt

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

