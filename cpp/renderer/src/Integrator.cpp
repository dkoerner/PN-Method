#include <integrator.h>











// returns sigma_t at sampled position (is invalid when we exceeded maxt)
double delta_tracking( const Scene* scene, const Ray3d& ray, double maxt, int component, RNGd& rng, V3d& sigma_t, bool debug )
{
	double sigma_t_max = scene->volume->getMaxExtinction()[component];

	double t = 0.0;
	while(true)
	{
		double step = -log( 1.0-rng.next1D() )/sigma_t_max;
		t += step;

		if(t>= maxt)
			break;

		sigma_t = scene->volume->evalExtinction(ray(t), debug);

		// russian roulette
		if(rng.next1D()<sigma_t[component]/sigma_t_max)
			break;
	}

	return t;
}






V3d nee( const TraceInfo& ti, RNGd& rng )
{
	// sample light source and transmittance ---
	LightSample ls;
	ls.refP = ti.current_vertex.getPosition();
	V3d attenuated_light_over_pdf = ti.scene->sample_attenuated_directlight( ls, rng );

	// apply scattering ---
	V3d phase_times_sigma_s = ti.scene->volume->evalPhase( ti.current_vertex.getPosition(), -ls.d, -ti.current_direction )*ti.current_vertex.m_sigma_s;

	return attenuated_light_over_pdf.cwiseProduct(phase_times_sigma_s).cwiseProduct(ti.throughput_over_pdf);
}





