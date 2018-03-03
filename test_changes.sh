docker run -it --mount type=bind,source=/Users/si992/code,target=/local d98112f6df0a
cp /local/guideseq/guideseq/guideseq.py /guideseq/guideseq/
cp /local/guideseq/guideseq/validation.py /guideseq/guideseq/
cp /local/guideseq/test/test_manifest.yaml /guideseq/test/
cp /local/guideseq/test/test_manifest_step2.yaml /guideseq/test/
rm -rf /guideseq/test/output/*
python /guideseq/guideseq/guideseq.py step1 -m /guideseq/test/test_manifest.yaml
python /guideseq/guideseq/guideseq.py step2 -m /guideseq/test/test_manifest_step2.yaml
