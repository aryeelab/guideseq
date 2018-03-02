docker run -it --mount type=bind,source=/Users/si992/code,target=/local d98112f6df0a
cp /local/guideseq/guideseq/guideseq.py /guideseq/guideseq/guideseq.py
cp /local/guideseq/test/test_manifest.yaml /guideseq/test/test_manifest.yaml
python /guideseq/guideseq/guideseq.py step1 -m /guideseq/test/test_manifest.yaml
