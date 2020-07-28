:: Script to convert presolved problems in .pre format to .mps
:: Author: Ivan Viviani

@echo off

setlocal enabledelayedexpansion
for %%f in (test_pre\*.pre) do ( 
	::echo "fullname: %%f"
    ::echo "name: %%~nf"
	cplex -c "read %%f" "write test_pre_mps\%%~nf_pre.mps"
)