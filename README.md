Parse the output from the Linux module avail command to generate an XSEDE compatible spreadsheet (CSV file) for updating the XSEDE software page

Usage: python gensoft.py > listing.csv

Notes:

1) I don’t currently do anything with /share/apps/compute/modulefiles. There’s no information on versions and there is some overlap with the output from the modules command. Figure that it will be easier to handle these by hand for now.

(2) You need to manually populate the python dictionaries with the type, domain, handle type and description. Let me know if you would like to just change the handle type to “module” for everything since that seems to be what you’re doing already.

(3) I execute a ‘module purge’ before running module avail so that we don’t pick up things like the MPI libraries and other Intel-dependent stuff.
