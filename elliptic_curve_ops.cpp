/*
 * powmod.cpp
 *
 *  Created on: Jun 7, 2018
 *      Author: jaydi
 */

#include <cstdint>
#include <iostream>
#include <cmath>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using namespace boost::multiprecision;

/*
 * ** Boost installation : https://www.boost.org/users/download/
 * ** just download, unzip, and use , no build required
 * ** cygwin : https://www3.ntu.edu.sg/home/ehchua/programming/howto/EclipseCpp_HowTo.html
 */



/*
	 * modulo_power calculates (base^exponent) % modulo , where exponent can be very large
	 * It exploits the following  formula
	 * (x * y)% z = ( x % z ) * ( y % z ), when x = y , it reduces to (x^2)%z = ( x % z ) * ( x % z ) =  ( x  % z )  ^2
	 * However ( x  % z ) ^2 may itself become more than z , so the final formula needs to be
	 * (x^2)%z = ( ( x  % z ) ^2  ) % z -- we refer this as the "repeated squaring formula"
	 * so the idea is express large exponent B as product of 2's powers and  carry  forward like below
	 * Example, base = 5, exponent = 117, modulo = 19,
	 *
	 * so ( 5 ^ 117 ) % 19
	 * = ( 5^64 % 19 )   *( 5^32 % 19 )  *( 5^16 % 19)    * ( 1 % 19)       * (   5^2 % 19 )   * ( 1 % 19 )   *   ( 5^1 % 19 )
	 *
	 * 5^64 % 19 = (( 5^32 ) % 19 )^2 % 19
	 * 5^32 % 19 = (( 5^16 ) % 19 )^2 % 19
	 * 5^16 % 19 = (( 5^8 ) % 19 )^2 % 19
	 * 5^8 % 19 = (( 5^4 ) % 19 )^2 % 19
	 * 5^4 % 19 = (( 5^2 ) % 19 )^2 % 19
	 * 5^2 % 19 = (( 5^1 ) % 19 )^2 % 19
	 *
	 * Exponent can be expressed as
	 * 		2^6*1		+	2^5*1		+	2^4*1			+	2^3*0			+	2^2*1		+	2^1*0	+		2^0*1 = 117
	 * 		 1        		1          		1            		0  					1				0				1	  : binary bits
	 *
	 * Positional value by	using the repeated squaring formula
	 * v =	(9*9)%19=5		(16*16)%19=9	(4*4)%19=16			(17*17)%19=4		(6*6)%19=17		(5*5)%19=6		5%19=5
	 *
	 * Accumulated Value : if the binary bit is on the accumulated value at that position=( positional value * prev accumulatd value)%modulo
	 * otherwise accumulated value does not change
	 * a = 	(5*4)%19=1		(9*11)%19=4		(16*9)%19=11		9					(17*5)%19=9		  5				(5*1) % 19 = 5
	 *
	 */
int512_t modulo_power(  int512_t base , int512_t exponent,  int512_t mod ) {
	int512_t  v , a=1;
	int counter = 0;
	while(exponent > 0 ) {
			v = ( counter ==0 ? base % mod : ((v % mod )* ( v % mod )) % mod );
			if ( exponent % 2 == 1 )
				a = (( a % mod ) *(  v % mod )) % mod;
			exponent /= 2;
			counter++;
	}
	return ( a< 0 ? a + mod : a );
}


/* calculates inverse of v such that ( inverse(v) * v ) % mod  = 1 */
/* Fermat's little theorem */
int512_t modulo_inverse(  int512_t v ,  int512_t mod ) {
	return( modulo_power( v, mod-2, mod ));
}


class el_point {
public:
	int512_t x,y;
	bool state ;
	el_point() {
		x=y=0;
		state=false;
	}
	el_point(int512_t _x, int512_t _y ) { x=_x; y=_y;state = true; }
	el_point( const el_point & rhs) { x = rhs.x ;y = rhs.y ; state = rhs.state;}
	el_point& operator=(const el_point& rhs)
	{
		if ( this != &rhs ) {x = rhs.x ; y = rhs.y ; state = rhs.state; }
		return *this ;
	}



	/**
	 * implements : this = 2*this
	 * y^2 = x^3 + ax + b
	 * c = (3*v.x^2 + a) / 2*v.y
	 * r.x = c^2 - 2*p.x
	 * r.y = c * (p.x - r.x) - p.y
	**/
	void doubler( int512_t a, int512_t mod ) {
		if ( !state ) return;
		int512_t num   = (3*( modulo_power(  x , 2,   mod )) + a ) % mod   ;
		if ( num < 0 ) { num += mod; }
		int512_t denom =  modulo_inverse( 2*y , mod );
		int512_t c     = (( num % mod )  * (denom % mod )) % mod;
		int512_t newx  = ( c *c  - 2 * x ) % mod;
		if (newx < 0 ) newx += mod;

		int512_t newy  = ( c * ( x - newx) - y )%mod;
		if (newy < 0 ) newy += mod;
		x = newx;
		y = newy;
	}

	/**
	 * implements this = this+q
	 * c = (q.y - p.y) / (q.x - p.x)
	 * r.x = c^2 - p.x - q.x
	 * r.y = c *(p.x - r.x) - p.y
	 **/
	void adder( const el_point & q, int512_t mod){
		if ( !q.state ) return;
		if ( !state ) { x = q.x; y = q.y ; state = true;
		return; }

		int512_t num   = (q.y - y) % mod;
		if ( num < 0 ) { num += mod; }

		int512_t denom = modulo_inverse( q.x - x, mod );
		if ( denom < 0 ) { denom += mod; }

		int512_t c = ((num %mod ) * (  denom % mod )) % mod;
		int512_t newx = ( c*c - x - q.x) % mod;
		if (newx < 0 ) newx += mod;

		int512_t newy = (c*(x - newx) -y ) % mod;
		if (newy < 0 ) newy += mod;

		x = newx;
		y = newy;
	}

	/* TODO squareroot */

} ;

/**
 * y^2 = x^3 + ax + b
 **/
class el_curve_ops
{
private :
	int512_t a,b, mod;
public :
	el_curve_ops( int512_t a_, int512_t b_, int512_t mod_ )
	{
		a= a_; b=b_; mod=mod_;
	}

	el_point scalarMult( int512_t privkey, const el_point & key)
	{
		if ( !key.state ) return key;
			el_point v(key), result;
			if ( privkey %2 == 1 )
			{ result=v; }
			privkey /= 2;
			while ( privkey > 0 ) {
				v.doubler(a,mod);
				if ( privkey %2 == 1)
					result.adder( v , mod);
				privkey /= 2;
			}
			return result;
	}


	bool verify(  const  el_point & v )
	{
		int512_t lhs= modulo_power(  v.y , 2,  mod );
		int512_t rhs1= modulo_power(  v.x , 3,  mod );
		int512_t rhs = (rhs1 + (a*v.x + b ))%mod ;
		return ( lhs == rhs );
	}
};

/**
 * modulo = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
** = 2^256 - 2^32 - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1
** G =
** ( 79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
     483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8 )
 */
int main()
{
  int512_t a = 0;
  int512_t b = 7;
  int512_t prime_modulo ( "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F");

  el_point generator_point( int512_t("0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798"),
		         int512_t("0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8" ));
  el_curve_ops eco( a,b,prime_modulo);

  //////////
  // matches with example from
  // https://bitcoin.stackexchange.com/questions/25024/how-do-you-get-a-bitcoin-public-key-from-a-private-key
  //////////
  //int512_t  privkey ("0x18E14A7B6A307F426A94F8114701E7C8E774E7F9A47E2C2035DB29A206321725");


  //////////
  // matches with example from
  // https://asecuritysite.com/encryption/bit_keys
  //////////
  //int512_t  privkey ("0xd8a8bb5aa721409deb930e8c2278b444d1bdb0f0a8a6e8cb97ec0ea9167175c5");

  // matches example from "mastering bitcoin book , page 78
  int512_t  privkey ("0x3aba4162c7251c891207b747840551a71939b0de081f85c4e44cf7c13e41daa6");

  cout << (eco.verify(generator_point )?" generator point is on elliptic curve":"generator point is NOT on elliptic curve") << endl;

  el_point pubkey = eco.scalarMult( privkey, generator_point );
  cout << "privkey="<< privkey<< " public key=> (" << pubkey.x << "," << pubkey.y  << ")" << (eco.verify(pubkey )?"pubkey point is on elliptic curve":"pubkey is NOT on elliptic curve") << endl;

}










