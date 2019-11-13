# Micro-QTL project

Scripts created for the micro-QTL project.


## create\_go\_tree.py

Converts the [GO flat file](http://purl.obolibrary.org/obo/go/go-basic.obo) to a pickled format optimised for tree traversal.

```sh
$ ./create_go_tree.py --output go_tree.pickle.xz go-basic.obo
```


## make\_uniprot\_idmapping\_db.py

Creates ID conversion maps from the [Uniprot idmapping file](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz).

```sh
# NCBI gene ID to GO-terms lookup table
$ ./make_uniprot_idmapping_db.py --from geneid --to go idmapping_selected.tab.gz output.pickle.xz
```
