=======
ContiGo
=======

SUMMARY
-------

ContiGo is a program that visualizes genome assemblies in a web browser. 

INSTALL
-------

ContiGo requires the Python Image Library (PIL). You can find it here:

http://www.pythonware.com/products/pil/

Download the latest version of ContiGo at:

https://github.com/vforget/Contigo/downloads

Extract the archive and install contiGo::

  $ tar xvzf Contigo-version.tar.gz
  $ cd Contigo-version
  $ python setup.py install

The last command may require root priviledges. As a result you may need to::

  $ sudo python setup.py install

USAGE
-----

NOTE: ContiGo works best with a recent version of Firefox. Google Chrome and Opera are functional, but buggy. Internet Explorer is not supported.

To generate a contiGo assembly view of an ACE file::

  $ contigo -i <path_to_ace> -o <out_dir>

For example::

  $ mkdir -p ~/tmp
  $ contigo -i ~/assembly_dir -o ~/tmp

ContiGo will seach ~/assembly_dir for a 454Contigs.ace file and put the output in ~/tmp.

To view the assembly, direct your web browser to ~/tmp/contigo.html.

To give it a try a small example ACE file is available here:

https://github.com/vforget/Contigo/tree/master/misc

By default contiGo searches for 454Contigs.ace. To use another ace file name use the -a option::

  $ contigo -i ~/assembly_dir -o ~/tmp -a some_other_assembly.ace

To test contiGo out on a larger assembly you can skip generating the pileup images with the -n option::

  $ contigo -i ~/assembly_dir -o ~/tmp -n

For more information on how to use contigo see the help message::

  $ contigo --help

EXAMPLES
--------

A live demo of contiGo for a genome assembly of *E. fergusonii* isolate ECD-227 can be found at:

http://www.genomequebec.mcgill.ca/compgen/contigo/ECD227/contigo.html

AUTHOR
------
Written by Vince Forgetta.
Contact: vincenzo.forgetta at mail.mcgill.ca
