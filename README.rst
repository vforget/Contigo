=======
ContiGo
=======

SUMMARY
-------

Contigo is a program that visualizes genome assemblies in a web browser. 


INSTALL
-------

ContiGo requires the Python Image Library (PIL). You can find it here::

  http://www.pythonware.com/products/pil/

To install Contigo, go to the directory where this README file is located and run::

  $ python setup.py install

This will require root priviledges. As a result you may need to::

  $ sudo python setup.py install

USAGE
-----

NOTE: Contigo works best with a recent version of Firefox. Google Chrome and Opera are functional, but buggy. Internet Explorer is not supported.

To generate a contiGo assembly view of an ACE file::

  python <contigo_path>/bin/contigo -i <path_to_ace> -o <out_dir>

For example, if you put contiGo in your home directory::

  $ mkdir -p ~/tmp
  $ python ~/Contigo/bin/contigo -i ~/Contigo/misc/ -o ~/tmp

Browse to ~/tmp/contigo.html in a web browser.

To test contiGo out on a larger assembly you can skip generating the pileup images with the -n option::

  $ python ~/Contigo/bin/contigo -i ~/Contigo/misc/ -o ~/tmp -n

By default contiGo searches for 454Contigs.ace. To use another ace file name use the -a option::

  $ python ~/Contigo/bin/contigo -i ~/Contigo/misc/ -o ~/tmp -a some_other_assembly.ace

For more information on how to use contigo see the help message::

  $ ~/Contigo/bin/contigo --help

EXAMPLES
--------

A live demo of contiGo for a genome assembly of *E. fergusonii* isolate ECD-227 can be found at::

  http://www.genomequebec.mcgill.ca/compgen/contigo/ECD227/contigo.html

AUTHOR
------
Written by Vince Forgetta.
Contact: vincenzo.forgetta at mail.mcgill.ca
