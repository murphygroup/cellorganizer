#!/bin/bash
# 
# Use module if available
#
# 2022-07-11
# Copyright 2022 Murphy Lab, CMU

function exit_code_to_boolean()
{
    [[ "$1" -ne "0" ]]
    echo "$?"
}

function module_if_available()
{
    which modulecmd >/dev/null 2>&1 ; module_available=$(exit_code_to_boolean)
    
    if (( module_available )); then
        module "$@"
    fi
}
