# waveletLib
Wavelet library for SDR#

This is a wavelet library plugin for SDR#, ported from the wavelet library developed by Ian Kaplan (www.bearcave.com).  It currently supports wavelet transforms based off the D4, Haar and Linear Interpolation basis functions.  It is meant to be used by application developers that need wavelet-based filters for their SDR# projects.

To use it, just copy the waveletLib directory and all of its sources into your SDR# SDRDEV subdirectory.  Then, include the library in your plugin source, and call the forwardTransform and reverseTransform functions for the appropriate basis.  

After compiling your plugin, copy both DLLs (the waveletLib DLL and the plugin you developed) into the SDRSHARP directory that houses your executable.  To register your plugin, edit the SDRSHARP XML file to add your plugin.  Now, your plugin should be loaded into SDR# the next time it starts.

