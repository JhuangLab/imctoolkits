alias bfg="java -jar ~/bin/tools/bfg-1.14.0.jar"
git add .
git commit -m "Lazy commit"
git commit --amend -CHEAD
bfg --strip-blobs-bigger-than 50M
git reflog expire --expire=now --all && git gc --prune=now --aggressive
git push
