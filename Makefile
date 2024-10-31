VERSION = "0.10.4"
PKG_NAME = 'zhunter'

upload2pip: pybuild twine

pybuild:
	@echo 'Building ${PKG_NAME}...'
	@python -m build

twine:
	@echo 'Uploading ${PKG_NAME}:${VERSION} to official pip'
	@python -m twine upload dist/*
	


