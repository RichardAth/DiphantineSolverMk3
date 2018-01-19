Solves diophantine equations of the form Ax² + Bxy + Cy² +Dx Ey +F = 0

This is based entirely on Dario Alpert's well-known solver.

His version was written as a web page using Java around 2003, and is still available. 
However modern browsers don't like Java. Firefox no longer support Java at all. I 
just could not get it to work with Internet Explorer, Edge or Chrome so I downloaded 
the source and converted it to a console program, which does work.

To make it work as a console program it was necessary to remove all the HTML tags, 
which means that the output is not as 'pretty' and the layout is a bit different.

Many other technical changes were necessary:

The code is very complex and difficult to understand.

In an effort to tidy it up I created a C++ equivalent program. In this program the 
functions were restructured much more radically and many comments added to make it 
a bit easier to follow. Also, I discovered that there were two completely different 
types of home-made bigintegers used, whch needed completely different functions to 
manipulate them. I replaced both types with GMP/MPIR extended precision functions, 
whch allowed me to simply remove a significant amount of code, and we can be confident 
that the extended-precision functions are well-documented and reliable. The division 
in particular was troublesome. It turned out that in some cases it is essential to 
use 'floor' division and in other cases 'truncation' division must be used.

I also grabbed Dario Alpert's web page describing his methods and converted it to a 
word-processor document.

I used Visual Studio, so non-windows users will need to adapt it.

This Mk3 version is functionally the same as the Mk2 version. It now uses the Boost 
Multi-precision library on top of the MPIR library. This allows the bigintegers to 
be handled more or less like normal integers. (The native MPIR functions are a bit 
like using assembler, Each operation would need a separate function call, whereas 
with Boost library normal expressions can be used, similar to Python)

A certain amount of tidying up has been done to make the code a bit less obscure 
e.g. reduce the usage of global variables to pass data around.

A useful feature of the Boost multi-precision library is that integers are 
automatically extended to multi-precision when necessary. Conversion from extended 
precision to normal integers requires a special function I wrote, which checks that 
the number will actually fit into 64 bits, otherwise the program aborts.

The method of factorising numbers has been changed to use a list of prime numbers,
and finds factors by trial division. This is much faster for large numbers, and
can now handle equations with 9 digit co-efficients in a reasonable amount of time.