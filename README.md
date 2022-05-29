# TDigest

This is a native Julia implementation of the t-digest that is at near parity with the Java
reference implementation. Some aspects of the code have been simplified relative to
the Java version, but the functionality should be essentially identical.

The current plan is to prepare this package for inclusion in the `OnlineStats` package.

# Things to do

In case you have a hankering to hack, please help with the following:

1. Port more tests, particularly those to do with duplicate points and accuracy testing
2. Add appropriate exports
3. Comment on the API design
4. Get some documentation started
5. Add CI/CD hooks using Github actions
