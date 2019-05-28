#!/bin/bash

cd `dirname "$0"`
cd "data"

test -e 'cancer' || mkdir cancer
cd cancer

for cancer in aml breast colon gbm kidney liver lung melanoma ovarian sarcoma
do
	test -e "$cancer" && continue
	mkdir "$cancer"
	cd "$cancer"
	wget "http://acgt.cs.tau.ac.il/multi_omic_benchmark/data/$cancer.zip"
	unzip "$cancer.zip"
	rm "$cancer.zip"
	cd ..
done

cd ..
test -e "clinical"  ||
	wget http://acgt.cs.tau.ac.il/multi_omic_benchmark/data/clinical.zip &&
	unzip clinical.zip &&
	rm clinical.zip

