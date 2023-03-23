VERSION = "0.10.4"
PKG_NAME = 'zhunter'

upload2pip: version pybuild twine

version:
	@echo 'Using version value defined in Makefile'
	@echo 'Updating to version ${VERSION} in following files:'
	@echo 'src/${PKG_NAME}/__init__.py'
	@./update_version.sh ${VERSION}

pybuild:
	@echo 'Building ${PKG_NAME}...'
	@python -m build

twine:
	@echo 'Uploading ${PKG_NAME}:${VERSION} to official pip'
	@python -m twine upload dist/*
	


