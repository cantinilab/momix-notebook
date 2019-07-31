#!/bin/bash

cd `dirname "$0"`
cd "data"

test -e 'bio_annotations' || mkdir bio_annotations
cd bio_annotations

MSIGBASE="http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2"

for annot in c2.cp.reactome.v6.2.symbols.gmt h.all.v6.2.symbols.gmt c5.all.v6.2.symbols.gmt
do
    test -e "$annot" || curl -o "$annot" "$MSIGBASE/$annot"
done

