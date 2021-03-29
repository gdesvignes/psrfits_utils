/* medfilt_new.c    Median filter in C*/
// algo by Phil Ekstrom
#ifndef EKSTROM_H
#define EKSTROM_H


#include <stdio.h>
#define STOPPER 0 /* Smaller than any datum */


template<class T> class Pair;
template<class T> class RunningMean;

template<class T> class Pair {

    public:
	Pair   *point;  /* Pointers forming list linked in sorted order */
	T value;	/* Values to sort */
};

template<class T> class RunningMean {
    private:
	unsigned int boxsize=10000;
	Pair<T> *buffer;
	Pair<T> *datpoint;  /* pointer into circular buffer of data */;
	Pair<T> small={NULL,STOPPER} ;  /* chain stopper. */
	Pair<T> big={&small,0} ;  /* pointer to head (largest) of linked list.*/

	Pair<T> *successor   ;  /* pointer to successor of replaced data item */
	Pair<T> *scan        ;  /* pointer used to scan down the sorted list */
	Pair<T> *scanold     ;  /* previous value of scan */
	Pair<T> *median;     ;  /* pointer to median */

    public:

	void set_boxsize(unsigned int bs) {boxsize = bs; buffer = new Pair<T>[boxsize](); datpoint=buffer;}

	T get_value(T datum) {

	    unsigned int i;

	    if(datum == STOPPER) datum = STOPPER + 1; /* No stoppers allowed. */
	    if( (++datpoint - buffer) >= boxsize) datpoint=buffer;  /* increment and wrap data in pointer.*/
	    datpoint->value=datum        ;  /* Copy in new datum */
	    successor=datpoint->point    ;  /* save pointer to old value's successor */
	    median = &big                ;  /* median initially to first in chain */
	    scanold = NULL               ;  /* scanold initially null. */
	    scan = &big                  ;  /* points to pointer to first (largest) datum in chain */
	    /* Handle chain-out of first item in chain as special case */
	    if( scan->point == datpoint ) scan->point = successor;
	    scanold = scan ;            /* Save this pointer and   */
	    scan = scan->point ;        /* step down chain */
	    /* loop through the chain, normal loop exit via break. */
	    for( i=0 ;i<boxsize ; i++ )
	    {
		/* Handle odd-numbered item in chain  */
		if( scan->point == datpoint ) scan->point = successor;  /* Chain out the old datum.*/
		if( (scan->value < datum) )        /* If datum is larger than scanned value,*/
		{
		    datpoint->point = scanold->point;          /* chain it in here.  */
		    scanold->point = datpoint;          /* mark it chained in. */
		    datum = STOPPER;
		};
	      /* Step median pointer down chain after doing odd-numbered element */
		median = median->point       ;       /* Step median pointer.  */
		if ( scan == &small ) break ;        /* Break at end of chain  */
		scanold = scan ;          /* Save this pointer and   */
		scan = scan->point ;            /* step down chain */
	      /* Handle even-numbered item in chain.  */
		if( scan->point == datpoint ) scan->point = successor; 
		if( (scan->value < datum) )         
		{
		    datpoint->point = scanold->point;       
		    scanold->point = datpoint;
		    datum = STOPPER;
		};
		if ( scan == &small ) break;
		scanold = scan ;                            
		scan = scan->point;
	    };
	    return( median->value );
	}
};

#endif
