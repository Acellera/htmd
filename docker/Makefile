DEST_DIR = '/app/'
current_dir = $(shell pwd)

build:
	cp ../package/htmd-deps/DEPENDENCIES .
	docker build -t htmdcontainer .

buildnocache:
	docker build -t htmdcontainer --no-cache=True .

run:
	docker run --mount type=bind,source=$(current_dir)/..,destination=$(DEST_DIR) --tty --interactive htmdcontainer

runtests:
	docker run -d --mount type=bind,source=$(current_dir)/..,destination=$(DEST_DIR) --tty --name htmdcontainer htmdcontainer
	docker container cp $(HOME)/Software/miniconda3/conda-bld/ htmdcontainer:/home/user/conda/
	docker exec -it htmdcontainer sh -c "conda install /home/user/conda/conda-bld/*-64/htmd-deps-[0-9]*.tar.bz2 -y -q"
	docker exec -it htmdcontainer sh -c "conda install /home/user/conda/conda-bld/*/htmd-data-[0-9]*.tar.bz2 -y -q"
	docker exec -it htmdcontainer sh -c "conda install /home/user/conda/conda-bld/*-64/htmd-[0-9]*.tar.bz2 -y -q"
	docker exec -it htmdcontainer sh -c "export HTMD_NONINTERACTIVE=1"
	docker exec -it htmdcontainer sh -c "python htmd/tests/run_inline_file_tests.py"
	docker exec -it htmdcontainer sh -c "htmd/tests/run_other_tests.sh"

stop:
	docker container stop htmdcontainer
	docker container rm htmdcontainer

rundetached:
	docker run -d htmdcontainer

push:
	docker tag htmdcontainer acellera/htmd:0.0.1
	docker push acellera/htmd:0.0.1

clean:
	docker image rm -f htmdcontainer
