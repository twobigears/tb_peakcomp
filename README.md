tb_peakcomp~
============

A stereo (L+R) peak compressor for Pd/libpd. Free = MIT licensed (use it, change it, etc).

This is an external for Pure Data (Pd), created to keep a check on the master stereo output. Useful for dynamic audio content especially on mobile devices.

Features:
* Threshold (dB)
* Ratio
* Attack (ms)
* Release (ms)
* Variable knee (hard knee to "smooth" knee)
* Makeup gain (dB)

Compiled for OSX and Android. Tested on OSX, Android and iOS (for iOS check the libpd for iOS Wiki).

Installation:
-------------

OSX: Add the contents of "OSX" to your Pd path
Android: Add the contents of the Android folder (all three folders) to the 'libs' folder of your Android project
iOS: Follow the steps in the libpd-for-iOS wiki

![tb_peakcomp~ demo patch](http://twobigears.com/othermedia/tb_peakcomp_demo_screen.png)


