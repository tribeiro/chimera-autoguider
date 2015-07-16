chimera-autoguider plugin
=======================

This is a telescope auto-guider plugin for the chimera observatory control system
(https://github.com/astroufsc/chimera).


This controller is reponsible for proving auto-guiding capabilities
to a telescope. It works by acquiring a star in a CCD and trying to
guarantee that it remain in its initial position, offsetting the
telescope accordingly. The operation is divided in two steps.

1) Finding a guide star
2) Guiding

Finding a guide star consists of:

1) Take an acquisition image
2) Check if a predetermined guide-star is in the FoV or try to
    find a suitable bright star for guiding. Raise an exception if
    no guide star is found.
3) Calculate the position of the star in the CCD (using centroid)

Guiding consists of

1) Take an image
2) Calculate the centroid of the guide-star
3) Offset the telescope
4) Repeat until stopped or an error occur.

Usage
-----



Installation
------------

Installation instructions.

::

   pip install -U chimera-autoguider

or

::

    pip install -U git+https://github.com/astroufsc/chimera_autoguider.git


Configuration Example
---------------------

Here goes an example of the configuration to be added on ``chimera.config`` file.

::

    instrument:
        name: model
        type: Example


Tested Hardware (for instruments)
---------------------------------

This plugin was tested on these hardware:

* Hardware example 1, model 2
* Hardware example 2, model 3


Contact
-------

For more information, contact us on chimera's discussion list:
https://groups.google.com/forum/#!forum/chimera-discuss

Bug reports and patches are welcome and can be sent over our GitHub page:
https://github.com/astroufsc/chimera_template/