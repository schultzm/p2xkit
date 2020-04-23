# p2xkit

`ispcr` Uses Biopython, emboss, primersearch to align primers to template to find amplimers (in-silico PCR).
`qpcr` does same as `ispcr` but adds qPCR probes to the mix (in-silico qPCR).

To use:  

```{bash}
p2xkit
usage: p2xkit [-h]  ...

optional arguments:
  -h, --help  show this help message and exit

Sub-commands help:
  
    ispcr     Perform in-silico PCR.
    qpcr      Perform in-silico qPCR.
    test      Run p2xkit unittests
    version   Print version to stdout
```

