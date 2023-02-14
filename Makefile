VERSION = "0.9.6"
version:
	@echo 'Using version value defined in Makefile'
	@echo 'Updating to version ${VERSION} in following files:'
	@echo 'src/zhunter/__init__.py'
	@./update_version.sh ${VERSION}