.PHONY: release

release:
ifndef version
	$(error usage: make release version=X.Y.Z)
endif
	git tag -fa "$(version)" -m "Release $(version)"
	git push origin "$(version)" --force
