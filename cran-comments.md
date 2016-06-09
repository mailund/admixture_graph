### Test environments
* i686-pc-linux-gnu, R 3.3.0 (Supposedly Educational)
* x86_64-w64-mingw32, R 3.2.2 (Fire Safety)

### R CMD check results
No ERRORS, WARNINGS or NOTES.
I have used the donttest command to skip most of the long examples during automated testing, but there's still three functions with examples taking more than 5 seconds (but less than 7). These three functions are pretty much the whole point of the package so having them tested sounds like the right thing to do, but I can skip them too if necessary.