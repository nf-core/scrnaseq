test: test_alevin test_star

test_alevin:
	nextflow run main.nf -profile test,docker

test_star:
	nextflow run main.nf -profile test,docker --aligner star

docker_build:
	docker build -t nfcore/scrnaseq:dev .

docker_push:
	sudo docker login
	sudo docker push nfcore/scrnaseq:dev
	docker images
