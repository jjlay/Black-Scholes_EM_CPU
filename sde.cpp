/*
 * JJ Lay
 * Middle Tennessee State University
 *
 * DATE        AUTHOR  COMMENTS
 * ----------  ------  ---------------
 * 2017-11-01  JJL     Initial version
 *
 * See:
 * 
 * Lay, H., Colgin, Z., Reshniak, V., Khaliq, A.Q.M. (2018).
 * "On the implementation of multilevel Monte Carlo simulation
 * of the stochastic volatility and interest rate model using
 * multi-GPU clusters." Monte Carlo Methods and Applications,
 * 24(4), pp. 309-321.
 * 
 * Lay, Harold A. (2020). "Simulation of stochastic systems using 
 * antithetic multilevel Monte Carlo on GPUs." Middle Tennessee
 * State University. Ph.D. Dissertation.
 * 
 */


//
// Standard Includes
//

#include <iostream>
#include <iomanip>
#include <random>


//
// Function: main()
//

int main(int argc, char *argv[]) {

	std::default_random_engine  gen;
	std::normal_distribution<double> dist(0.0, 1.0);


	// The simulation is executed for an increasingly greater number of 
	// steps (thus the step size becomes smaller).

	for (unsigned int steps = 2; steps < 10000; steps = steps * 2) {
		// METASAMPLES is the number of solutions used to estimate the
		// high-level statistics
		const int metasamples = 2000;

		// SAMPLES is the number of samples used to estimate the solution
		const unsigned int samples = 2000;

		// T is the time to expiry
		const double T = 1.0;

		// X0 is the initial value of the asset
		const double X0 = 100.0;

		// R is the constant interest rate
		const double r = 0.05;

		// VOLATILITY is the constanrt volatility of the asset
		const double volatility = 0.03;
		
		// Analytical solution
		const double analytical = X0 * exp(r * T);

		// Step size
		double dt = T / static_cast<double>(steps);
		double sqrtdt = sqrt(dt);

		// DATAEM captures the individual samples within a simulation
		auto dataEM = new double[samples];
		auto metameanEM = new double[metasamples];
		auto metastdevEM = new double[metasamples];

		auto dataMilstein = new double[samples];
		auto metameanMilstein = new double[metasamples];
		auto metastdevMilstein = new double[metasamples];

		double XEM = 0.0, XMilstein = 0.0, dx = 0.0, dW = 0.0;

		for (auto m = 0; m < metasamples; m++) {
			//
			// Perform simulation
			//
		
			for (auto i = 0; i < samples; i++) {
				XEM = X0;
				XMilstein = X0;

				for (auto j = 0; j < steps; j++) {
					// Use the same random values for both the Milstein and Euler-Maruyama
					dW = dist(gen) * sqrtdt;

					// Euler-Maruyama
					dx = (r * XEM * dt) + (volatility * XEM * dW);
					XEM += dx;

					// Milstein
					dx = (r * XMilstein * dt) + (volatility * XMilstein * dW) + (0.5 * (volatility * XMilstein) * (volatility) * (dW * dW - dt));
					XMilstein += dx;
				}
	
				dataEM[i] = XEM;
				dataMilstein[i] = XMilstein;
			}

			//
			// Calculate statistics
			//
	
			double meanEM = 0.0, meanMilstein = 0.0;
	
			for (auto i = 0; i < samples; i++) {
				meanEM += dataEM[i];
				meanMilstein += dataMilstein[i];
			}
	
			meanEM = meanEM / static_cast<double>(samples);
			meanMilstein = meanMilstein / static_cast<double>(samples);

			double stdevEM = 0.0, stdevMilstein = 0.0;

			for (auto i = 0; i < samples; i++) {
				stdevEM += pow(dataEM[i] - meanEM, 2.0);
				stdevMilstein += pow(dataMilstein[i] - meanMilstein, 2.0);
			}

			stdevEM = sqrt(stdevEM / static_cast<double>(samples));
			stdevMilstein = sqrt(stdevMilstein / static_cast<double>(samples));

			metameanEM[m] = meanEM;
			metastdevEM[m] = stdevEM;
			metameanMilstein[m] = meanMilstein;
			metastdevMilstein[m] = stdevMilstein;
		}

		double meanmean = 0.0, meanstdev = 0.0, stdevmean = 0.0, stdevstdev = 0.0;

		//
		// Euler-Maruyama
		//
		
		for (auto i = 0; i < metasamples; i++) {
			meanmean += metameanEM[i];
			stdevmean += metastdevEM[i];
		}

		meanmean = meanmean / static_cast<double>(metasamples);
		stdevmean = stdevmean / static_cast<double>(metasamples);

		for (auto i = 0; i < metasamples; i++) {
			meanstdev += pow(metameanEM[i] - meanmean, 2.0);
			stdevstdev += pow(metastdevEM[i] - stdevmean, 2.0);
		}

		meanstdev = sqrt(meanstdev / static_cast<double>(metasamples));
		stdevstdev = sqrt(stdevstdev / static_cast<double>(metasamples));

		std::cout << "Euler-Maruyama : dt : " << std::fixed << std::setw(6) << std::setprecision(5) << dt << " | "
			<< "MEAN :: mean: " << std::fixed << std::setw(8) << std::setprecision(7) << meanmean << " : stdev: " << meanstdev << " | "
			<< "STD DEV :: mean: " << stdevmean << " : stdev: " << stdevstdev << std::endl;
		
		//
		// Milstein
		//

		meanmean = 0.0;
		meanstdev = 0.0;
		stdevmean = 0.0;
		stdevstdev = 0.0;

		for (auto i = 0; i < metasamples; i++) {
			meanmean += metameanMilstein[i];
			stdevmean += metastdevMilstein[i];
		}

		meanmean = meanmean / static_cast<double>(metasamples);
		stdevmean = stdevmean / static_cast<double>(metasamples);

		for (auto i = 0; i < metasamples; i++) {
			meanstdev += pow(metameanMilstein[i] - meanmean, 2.0);
			stdevstdev += pow(metastdevMilstein[i] - stdevmean, 2.0);
		}

		meanstdev = sqrt(meanstdev / static_cast<double>(metasamples));
		stdevstdev = sqrt(stdevstdev / static_cast<double>(metasamples));

		std::cout << "Milstein       : dt : " << std::fixed << std::setw(6) << std::setprecision(5) << dt << " | "
			<< "MEAN :: mean: " << std::fixed << std::setw(8) << std::setprecision(7) << meanmean 
            << " : stdev: " << meanstdev 
            << " | "
			<< "STD DEV :: mean: " << stdevmean 
            << " : stdev: " << stdevstdev << std::endl;
	}
		
	return 0;
}

