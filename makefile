
main: 	README.md style.css
	pandoc README.md -f markdown -t html -o README.html --toc --from markdown-yaml_metadata_block --css=style.css --highlight-style=tango
