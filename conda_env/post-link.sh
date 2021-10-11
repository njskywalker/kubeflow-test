#!/usr/bin/env bash

if [ "$(uname)" == "Darwin" ]; then
    echo "Installing gfortran for macOS"
    conda install -c brittainhard gfortran
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    echo "Installing gfortran for Linux"
    conda install -c anaconda gfortran_linux-64
fi

jupyter-nbextension enable --py --user widgetsnbextension
jupyter-nbextension enable --py --user nglview