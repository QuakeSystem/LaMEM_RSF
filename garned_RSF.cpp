bool SECOND_ORDER_BCS_D = true;
bool SECOND_ORDER_BCS_I = true;
#include "garnet.h"

/*
STEP 9: Lapusta's algorithm
Step 9a: Include adaptive grids to BP3 from benchmark.cpp
Step 9d: Use two-step solver inspired by Lapusta 2000
       , use JÎµ = Î²_ * J * (Ï„-Ï„0) / Î¼ instead of Ï„
Step 9e: Put v at 0b11
Step 9h: remove iteration during coseismic
      ,  output only from rank0 failed
			,  use 40km*40km domain
			,  use same axis configuration for x and z
Renamed as rsf_2d_fd
// Add Rayleigh damping - removed
Rewrite aging law
Add vratio to control timestep
All velocity-weakening
*/

struct vNode {
	static const unsigned long value = 0b11;
};

static const size_t N = 2;

using namespace units;

const double Wf     =  20    *  km;
const double H      =  15    *  km;
const double h      =  3     *  km;
const double Ï      =  2670  *  kg/m3;
const double cs     =  3.464 *  km/s;
const double Î¼      =  Ï*cs*cs; // =~ 32e9

const double H0     =  /*2 */  Wf;
const double H1     =  /*2 */  Wf;
const double Hf     =  H + h/*2*/;
const double Hg     =  1 * km;

const double Vp     =  1e-9  *  m/s; // loading rate
const double V0     =  1e-6  *  m/s; // reference slip rate
const double Ïƒn     =  50    *  MPa;
const double Dc     =  0.008 *  m  ;

const double a0     =  0.010;
const double a_max  =  0.025;
const double b0     =  0.015;
const double f0     =  0.6  ;
// const double Î·      =  Î¼ / cs / 2.; // damping parameter

const double Ï„0     =  Ïƒn * a_max * std::asinh( Vp / (2*V0) * std::exp( ( f0 + b0 * std::log( V0 / Vp ) ) / a_max ) );

template< class B, class X >
void ComputeHorizontalAxis( Domain<N>& Î©, long int n1, long int nf, B& Î², X& x )
{
	double A0 = ( n1*n1*n1*n1 * ( H0 * nf - Hg * n1 ) ) / ( std::pow( n1 - nf, 4 ) * nf );
	double A1 = ( 4 * n1*n1*n1 * ( Hg * n1 - H0 * nf ) ) / std::pow( n1 - nf, 4 );
	double A2 = ( 6 * n1*n1 * nf * ( H0 * nf - Hg * n1) ) / std::pow( n1 - nf, 4 );
	double A3 = ( n1 * ( -4 * H0 * nf*nf*nf*nf + Hg * ( n1*n1*n1*n1 - 4 * n1*n1*n1 * nf + 6 * n1*n1 * nf*nf + nf*nf*nf*nf ) ) ) / ( std::pow( n1 - nf, 4 ) * nf );
	double A4 = ( nf*nf*nf * ( H0 * nf - Hg * n1 ) ) / std::pow( n1 - nf, 4 );

	// std::cout << "A0 = " << A0 << ", A1 = " << A1 <<  ", A2 = " << A2 <<  ", A3 = " << A3 <<  ", A4 = " << A4 <<  std::endl;

	Domain<N-1> Î“x   =  Î©.Face(-2);
	x.Set( [&]( Coor<N-1> Î¾ ) {
		double Î¶ = /*1 -*/ Î¾[0];
		return ( Î¶ <= double(nf) / double(n1) ?
			(n1*Hg)/nf * Î¶ :
			A0*Î¶*Î¶*Î¶*Î¶ + A1*Î¶*Î¶*Î¶ + A2*Î¶*Î¶ + A3*Î¶ + A4 );//
	}, Î“x.X() );
	x.ViewComplete();

	Î²[0][0].Set( [&]( Coord<N> Î¾ ) { // m
		double Î¶ = /*1 -*/ Î¾[0];
		return Î¶ <= double(nf) / double(n1) ?
			n1 * Hg / nf :
			4*A0*Î¶*Î¶*Î¶ + 3*A1*Î¶*Î¶ + 2*A2*Î¶ + A3;
	}, Î©.X() );
	Î²[0][0].template SetBC<Neumann>( -1, 0. );
	Î²[0][0].node(0b10).prescribe_bc(-1);
	Î²[0][0].node(0b11).prescribe_bc(-1);
	Î²[0][0].template SetBC<Neumann>( +1, 0. );
	Î²[0][0].node(0b10).prescribe_bc(+1);
	Î²[0][0].node(0b11).prescribe_bc(+1);
	Î²[0][0].template SetBC<Neumann>( -2, 0. );
	Î²[0][0].node(0b10).prescribe_bc(-2);
	Î²[0][0].node(0b11).prescribe_bc(-2);
	Î²[0][0].template SetBC<Neumann>( +2, 0. );
	Î²[0][0].node(0b10).prescribe_bc(+2);
	Î²[0][0].node(0b11).prescribe_bc(+2);
	Î²[0][0].ViewComplete();
}


template< class B, class Z, class DZ >
void ComputeVerticalAxis( Domain<N>& Î©, long int n1, long int nf, B& Î², Z& z, DZ& dz )
{
	double A0 = ( n1*n1*n1*n1 * ( H1 * nf - Hf * n1 ) ) / ( std::pow( n1 - nf, 4 ) * nf );
	double A1 = ( 4 * n1*n1*n1 * ( Hf * n1 - H1 * nf ) ) / std::pow( n1 - nf, 4 );
	double A2 = ( 6 * n1*n1 * nf * ( H1 * nf - Hf * n1) ) / std::pow( n1 - nf, 4 );
	double A3 = ( n1 * ( -4 * H1 * nf*nf*nf*nf + Hf * ( n1*n1*n1*n1 - 4 * n1*n1*n1 * nf + 6 * n1*n1 * nf*nf + nf*nf*nf*nf ) ) ) / ( std::pow( n1 - nf, 4 ) * nf );
	double A4 = ( nf*nf*nf * ( H1 * nf - Hf * n1 ) ) / std::pow( n1 - nf, 4 );

	//std::cout << "A0 = " << A0 << ", A1 = " << A1 <<  ", A2 = " << A2 <<  ", A3 = " << A3 <<  ", A4 = " << A4 <<  std::endl;

	Domain<N-1> Î“z   =  Î©.Face(-1);
	z.Set( [&]( Coor<N-1> Î¾ ) {
		double Î¶ = /*1 -*/ Î¾[0];
		return ( Î¶ <= double(nf) / double(n1) ?
			(n1*Hf)/nf * Î¶ :
			A0*Î¶*Î¶*Î¶*Î¶ + A1*Î¶*Î¶*Î¶ + A2*Î¶*Î¶ + A3*Î¶ + A4 );//
	}, Î“z.X() );
	z.ViewComplete();

	dz = z.D(0) / n1;

	Î²[1][1].Set( [&]( Coord<N> Î¾ ) { // m
		double Î¶ = /*1 -*/ Î¾[1];
		return Î¶ <= double(nf) / double(n1) ?
			n1 * Hf / nf :
			4*A0*Î¶*Î¶*Î¶ + 3*A1*Î¶*Î¶ + 2*A2*Î¶ + A3;
	}, Î©.X() );
	Î²[1][1].template SetBC<Neumann>( -1, 0. );
	Î²[1][1].node(0b10).prescribe_bc(-1);
	Î²[1][1].node(0b11).prescribe_bc(-1);
	Î²[1][1].template SetBC<Neumann>( +1, 0. );
	Î²[1][1].node(0b10).prescribe_bc(+1);
	Î²[1][1].node(0b11).prescribe_bc(+1);
	Î²[1][1].template SetBC<Neumann>( -2, 0. );
	Î²[1][1].node(0b10).prescribe_bc(-2);
	Î²[1][1].node(0b11).prescribe_bc(-2);
	Î²[1][1].template SetBC<Neumann>( +2, 0. );
	Î²[1][1].node(0b10).prescribe_bc(+2);
	Î²[1][1].node(0b11).prescribe_bc(+2);
	Î²[1][1].ViewComplete();
}


void RateAndState()
{
	////////////////////////////
	//                        //
	//   DOMAIN DEFINITION    //
	//                        //
	////////////////////////////

	Domain<N>   Î©( "domain" );
	Î©.Recenter( 0.5, 0.5 );

	int n0   =  GetOption<int>( "-domain_da_grid_x" ) - 1;
	int ng   =  GetOption<int>( "-domain_da_grid_g" ) - 1;
	int n1   =  GetOption<int>( "-domain_da_grid_y" ) - 1;
	int nf   =  GetOption<int>( "-domain_da_grid_f" ) - 1;

	if( nf > n1 ) throw std::logic_error( "nf > n1!" );


	////////////////////////////
	//                        //
	//   FIELDS               //
	//                        //
	////////////////////////////

	auto   v    =   NewODE<BDF<2>,AM<2>,AB<2>>( Î©.NewScalar(  "v",   vNode() ) );
	auto  JÎµ    =   NewODE<BDF<2>,AM<2>,AB<2>>(     Storage( "JÎµ",Grad(v(0)) ) );
	auto   Î¸    =   NewODE<BDF<2>,AM<2>,AB<2>>(     Storage(  "Î¸",     v(0).Face(-1) ) );
	auto   Ïƒ    =                                   Storage(  "Ïƒ",    JÎµ(0) );

	auto   U    =                                   Storage(  "U",   Î¸(0) );
	auto   V    =   NewODE<BDF<2>,AM<2>,AB<2>>(     Storage(  "V",   Î¸(0) ) );
	auto   Î¤    =                                   Storage(  "Î¤",   Î¸(0) );

	auto   Ï‡    =    NewTimekeeper( Î¸, v, JÎµ, V );


	////////////////////////////
	//                        //
	//   GEOMETRY             //
	//                        //
	////////////////////////////

	auto   Î²    =   Storage(  "Î²",   Prod<2>( Interpolate(JÎµ(0)), Interpolate(JÎµ(0)) ), Diagonal2ndOrder() );
	auto   g_   =   Storage(  "g_",  Prod<2>( Interpolate(JÎµ(0)), Interpolate(JÎµ(0)) ), Diagonal2ndOrder() );
	auto   J    =   Storage(  "J",   Î² * Î² );
	auto   x    =   Storage(  "x",   Î²[0][0].Face(-2) );
	auto   z    =   Storage(  "z",   Î²[1][1].Face(-1) );
	auto  dz    =   Storage( "dz",   Î¸(0) );

	ComputeHorizontalAxis( Î©, n0, ng, Î², x );//
	ComputeVerticalAxis  ( Î©, n1, nf, Î², z, dz );

	// g = Prod<2>( T(Î²), Î² );  // <--- segfault, check!!
	g_[0][0] = 1. / Î²[0][0] / Î²[0][0];
	g_[1][1] = 1. / Î²[1][1] / Î²[1][1];

	J = 1. / Sqrt( Det( g_ ) );
	J.SetBC<Neumann>( 0. );//


	////////////////////////////
	//                        //
	//   FAULT PARAMETERS     //
	//                        //
	////////////////////////////

	auto   a    =   Storage(  "a",   Î¸(0) );
	auto   b    =   Storage(  "b",   Î¸(0) );
	auto   Î¾    =   Storage(  "Î¾",   Î¸(0) ); // Lapusta et al. (2000, JGR)

	// a. Set( [&]( double z ) {
	// 	return std::min( a_max, std::max( a0, a0 + ( z - H ) * ( a_max - a0 ) / h ) );
	// }, z );
	a = a0;

	b = b0;

	// Lapusta et al. (2000, JGR)
	Î¾. Set( [&]( double a, double b, double dz )
	{
		double  A  =  a * Ïƒn;
		double  B  =  b * Ïƒn;
		double  Î³  =  pi / 4;
		double  k  =  Î³ * Î¼ / dz;
		double  Ï‡  =  0.25 * std::pow( k*Dc/A - (B-A)/A, 2 ) - k*Dc/A;
		double  Î¾1 =  std::min(  A / ( k*Dc   - (B-A) ), 0.5 );
		double  Î¾2 =  std::min( 1 - (B-A)/k/Dc, 0.5 );
		return  Ï‡ > 0 ? Î¾1 : Î¾2;
	}, a, b, dz );

	auto timestep_factor = GetOption<double>( "-timestep_factor" );
	Î¾ = Î¾ * timestep_factor;
	Î¾.ViewComplete();


	////////////////////////////
	//                        //
	//   OUTPUTTING           //
	//                        //
	////////////////////////////

	int i = 0;

	auto stdout_line = [&]() {
		std::cout <<    "i  = "  << i
		          <<  ", t  = "  << Ï‡.t() / yr << " yr"
		          <<  ", dt = "  << Ï‡.dt << " s"
		          //<<  ", U  = "  << Max(U) << " m"
		          <<  ", V  = (" << Min(V(0)) << ", " << Max(V(0)) << ")"
		          <<  ", Î¸  = (" << Min(Î¸(0)) << ", " << Max(Î¸(0)) << ")"
		          <<  ", Î¤  = (" << Min(Î¤) << ", " << Max(Î¤) << ")"
		          << std::endl;
	};

	auto json_line = [&]() {
		return json {
			{ "i",  i },
			{ "t",  Ï‡.t() },
			{ "dt", Ï‡.dt },
			{ "Vmax", Max(V(0)) },
			{ "Vmin", Min(V(0)) },
			{ "Î¤max", Max(Î¤) },
			{ "Î¤min", Min(Î¤) }
		};
	};

	std::fstream outfile;

	auto json_file_start = [&]() {
		outfile = std::fstream( "./output.json", std::ios_base::out | std::ios_base::trunc );
		//outfile << "{\"time_series\":[\n";
		outfile   << "\t" << json_line() << "\n";
		outfile.flush();
	};

	auto json_file_open = [&]() {
		outfile = std::fstream( "./output.json", std::ios_base::in | std::ios_base::out | std::ios_base::ate );
	};

	auto json_line_add = [&]() {
		//outfile.seekp( -4, std::ios_base::cur );
		outfile   << "\t" << json_line() << "\n";
		outfile.flush();
	};

	auto output_fields_1 = [&]() {
		Î¤.ViewComplete(i);
		V(0).ViewComplete(i);
		U.ViewComplete(i);
		Î¸(0).ViewComplete(i);
	};

	auto output_fields_2 = [&]() {
		//JÎµ(0).ViewComplete(i);
		Ïƒ.ViewComplete(i);
		v(0).ViewComplete(i);
	};


	////////////////////////////
	//                        //
	//   INITIAL CONDITIONS   //
	//                        //
	////////////////////////////

	auto  cp    =  Î©.NewCheckpointer( Name(Ï‡), Name(i) );

	v(0)	=   Vp / 2.;
	v(-1)	=   Vp / 2.;

	V(0)  =   Vp;
	V(-1) =   Vp;
	U     =   0;

	JÎµ(0) =   0;
	JÎµ(-1)=   0;
	Ïƒ     =   0;

	Î¸(0)  =   (Dc/V0) * Exp( (a/b) * Log( 2 * V0 / Vp * Sinh( Ï„0 / (a*Ïƒn) ) ) - f0 / b );
	Î¤     =   Ï„0;

	// Lapusta et al. (2000, JGR)
	Ï‡.dt = AbsMin( Î¾ * Dc / V(0) );
	Ï‡.dt = 100 * s;

	if( CheckOption( "-restart_checkpoint" ) ) {
		cp.Load( GetOption<int>( "-restart_checkpoint" ) );
		json_file_open();
	} else
		json_file_start();

	output_fields_1();
	output_fields_2();
	stdout_line();

	Ï‡.Step();

	i++;


	////////////////////////////
	//                        //
	//   PHYSICAL EQUATIONS   //
	//                        //
	////////////////////////////

	JÎµ.SetBC<Neumann >( 0. );

	v.SetBC<Dirichlet>( +2,  Vp/2 );
	v.SetBC<Dirichlet>( -1,  If( z < Wf, V(0)/2, Vp/2 ) );
	v.SetBC<Neumann  >( -2,  0. );
	v.SetBC<Neumann  >( +1,  0. );

	V.SetRHS(0);

	auto S = NewSolver( "benchmark", v(0) );
	// auto damp = GetOption<double>( "-damp_b" );

	S.SetResidual( [&]( auto& Rv ) {
		JÎµ.TrivialSolve<BDF<2>>( 0., g_ * J * Grad( v(0) ) );
		// Rv = v.template Residual<BDF<2>>( Div( JÎµ(0) ) / J * Î¼ / Ï );
		Rv = v.template Residual<BDF<2>>( Div( JÎµ(0) /*+ damp * Ï‡.dt * JÎµ.Differentiate<BDF<2>>()*/ ) / J * Î¼ / Ï );
	} );


	////////////////////////////
	//                        //
	//   MAIN LOOP            //
	//                        //
	////////////////////////////

	auto nsteps = GetOption<int>( "-n_steps" );
	auto output_interval_1 = GetOption<int>( "-output_interval_1" );
	auto output_interval_2 = GetOption<int>( "-output_interval_2" );
	auto checkpoint_interval = GetOption<int>( "-checkpoint_interval" );

	auto C = GetOption<double>( "-cfl_c" );

	for( ; i <= nsteps; )
	{
		v(0).Set(v(-1));
		V(0).Set(V(-1));

		// Ï‡.dt = AbsMin( Î¾ * Dc / V(0) );

		// Î¸.TrivialSolve<BDF<2>>( -V(0) / Dc , 1 );
		Î¸(0) = Dc / V(0) * ( 1 - ( 1 - V(0) / Dc * Î¸(-1) ) * Exp( -V(0) / Dc * Ï‡.dt ) );
		v(0).SetBC<Dirichlet>( -1,  If( z < Wf, V(0)/2, Vp/2 ) );

		if( cs * Ï‡.dt / (Hf/nf) < C ) {
			JÎµ.TrivialSolve<BDF<2>>( 0., g_ * J * Grad( v(0) ) );
			v .TrivialSolve<BDF<2>>( 0., Div( JÎµ(0) ) / J * Î¼ / Ï );
		} else {
			//S.Solve();
		}

		std::cout << "Step 1: " << "MinMax(v) = " << Min(v(0))  << ", " << Max(v(0)) << std::endl;
		std::cout << "        " << "MinMax(V) = " << Min(V(0))  << ", " << Max(V(0)) << std::endl;

		Ïƒ =  Î¼ * Î² * JÎµ(0) / J;

		V(0) =  If( z < Wf, ( V(-1) + Sinh( ( Ïƒ[0].Face(-1) + Ï„0 ) / (a * Ïƒn) ) * (2*V0) / Exp( ( f0 + b * Log( V0 * Î¸(0) / Dc ) ) / a ) )/ 2, Vp );

		// Î¸.TrivialSolve<BDF<2>>( -V(0) / Dc , 1 );
		Î¸(0) = Dc / V(0) * ( 1 - ( 1 - V(0) / Dc * Î¸(-1) ) * Exp( -V(0) / Dc * Ï‡.dt ) );
		v(0).SetBC<Dirichlet>( -1,  If( z < Wf, V(0)/2, Vp/2 ) );

		if( cs * Ï‡.dt / (Hf/nf) < C ) {
			JÎµ.TrivialSolve<BDF<2>>( 0., g_ * J * Grad( v(0) ) );
			v .TrivialSolve<BDF<2>>( 0., Div( JÎµ(0) ) / J * Î¼ / Ï );
		} else {
			JÎµ.TrivialSolve<BDF<2>>( 0., g_ * J * Grad( v(0) ) );
			v .TrivialSolve<BDF<2>>( 0., Div( JÎµ(0) ) / J * Î¼ / Ï );
			S.Solve();
		}

		// if( !S.Success() ) {
		// 	std::cout << "Unsuccessful solve, exiting!" << std::endl;
		// 	break;
		// }

		std::cout << "Step 2: " << "MinMax(v) = " << Min(v(0))  << ", " << Max(v(0)) << std::endl;
		std::cout << "        " << "MinMax(V) = " << Min(V(0))  << ", " << Max(V(0)) << std::endl;

		U += Ï‡.dt * V(0);

		Ïƒ =  Î¼ * Î² * JÎµ(0) / J;
		V(0) =  If( z < Wf, Sinh( ( Ïƒ[0].Face(-1) + Ï„0 ) / (a * Ïƒn) ) * (2*V0) / Exp( ( f0 + b * Log( V0 * Î¸(0) / Dc ) ) / a ), Vp );
		Î¤ = Ïƒ[0].Face(-1) + Ï„0;

		// double vratio = std::abs( Max(V(0)) / Max(V(-1)) );
		// if( vratio > 1e1 || vratio < 1e-1 ) {
		// 	std::cout << "V(0)/V(-1) = " << vratio << ", " << "Half time step and try again!" << std::endl;
		// 	Ï‡.dt = Ï‡.dt / 2;
		// 	continue;
		// }

		double dt_temp = Ï‡.dt;
		Ï‡.dt = AbsMin( Î¾ * Dc / V(0) );
		Ï‡.dt = std::min( Ï‡.dt, 1.2 * dt_temp );

		json_line_add();
		stdout_line();
		if( !( i % output_interval_1 ) )
			output_fields_1();
		if( !( i % output_interval_2 ) )
			output_fields_2();
		if( !( i % checkpoint_interval ) )
			cp.Save(i);

		Ï‡.Step();
		i++;
	}
}

int main( int argc, char **argv )
{
	InitGARNET( argc, argv );
	RateAndState();
	ExitGARNET();
}