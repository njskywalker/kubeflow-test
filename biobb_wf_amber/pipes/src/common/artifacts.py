from kfp import dsl

class PDBFile(dsl.Artifact):
    TYPE_NAME = "PDBFile"
    EXTENSION = ".pdb"
