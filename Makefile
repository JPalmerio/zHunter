VERSION = "0.10.2"

upload2pip: version pybuild twine

version:
	@echo 'Using version value defined in Makefile'
	@echo 'Updating to version ${VERSION} in following files:'
	@echo 'src/zhunter/__init__.py'
	@./update_version.sh ${VERSION}

pybuild:
	@echo 'Building zhunter...'
	@python -m build

twine:
	@echo 'Uploading to official pip'
	@python -m twine upload dist/*
	


