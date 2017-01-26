#pragma once

#include <vector>


/**
 * \brief Discrete probability distribution
 *
 * This data structure can be used to transform uniformly distributed
 * samples to a stored discrete probability distribution.
 *
 * \ingroup libcore
 */
template<typename T>
struct DiscretePDF
{
public:
	/// Allocate memory for a distribution with the given number of entries
	explicit inline DiscretePDF(size_t nEntries = 0)
	{
		reserve(nEntries);
		clear();
	}

	/// Clear all entries
	inline void clear() {
		m_cdf.clear();
		m_cdf.push_back(T(0.0));
		m_normalized = false;
	}

	/// Reserve memory for a certain number of entries
	inline void reserve(size_t nEntries) {
		m_cdf.reserve(nEntries+1);
	}

	/// Append an entry with the specified discrete probability
	inline void append(T pdfValue) {
		m_cdf.push_back(m_cdf[m_cdf.size()-1] + pdfValue);
	}

	/// Return the number of entries so far
	inline size_t size() const {
		return m_cdf.size()-1;
	}

	/// Access an entry by its index
	inline T operator[](size_t entry) const {
		return m_cdf[entry+1] - m_cdf[entry];
	}

	/// Have the probability densities been normalized?
	inline bool isNormalized() const {
		return m_normalized;
	}

	/**
	 * \brief Return the original (unnormalized) sum of all PDF entries
	 *
	 * This assumes that \ref normalize() has previously been called
	 */
	inline T getSum() const {
		return m_sum;
	}

	/**
	 * \brief Return the normalization factor (i.e. the inverse of \ref getSum())
	 *
	 * This assumes that \ref normalize() has previously been called
	 */
	inline T getNormalization() const {
		return m_normalization;
	}

	/**
	 * \brief Normalize the distribution
	 *
	 * \return Sum of the (previously unnormalized) entries
	 */
	inline T normalize() {
		m_sum = m_cdf[m_cdf.size()-1];
		if (m_sum > 0) {
			m_normalization = T(1.0) / m_sum;
			for (size_t i=1; i<m_cdf.size(); ++i)
				m_cdf[i] *= m_normalization;
			m_cdf[m_cdf.size()-1] = 1.0f;
			m_normalized = true;
		} else {
			m_normalization = T(0.0);
		}
		return m_sum;
	}

	/**
	 * \brief %Transform a uniformly distributed sample to the stored distribution
	 *
	 * \param[in] sampleValue
	 *     An uniformly distributed sample on [0,1]
	 * \return
	 *     The discrete index associated with the sample
	 */
	inline size_t sample(T sampleValue) const {
		typename std::vector<T>::const_iterator entry =
				std::lower_bound(m_cdf.begin(), m_cdf.end(), sampleValue);
		size_t index = (size_t) std::max((ptrdiff_t) 0, entry - m_cdf.begin() - 1);
		return std::min(index, m_cdf.size()-2);
	}

	/**
	 * \brief %Transform a uniformly distributed sample to the stored distribution
	 *
	 * \param[in] sampleValue
	 *     An uniformly distributed sample on [0,1]
	 * \param[out] pdf
	 *     Probability value of the sample
	 * \return
	 *     The discrete index associated with the sample
	 */
	inline size_t sample(T sampleValue, T &pdf) const {
		size_t index = sample(sampleValue);
		pdf = operator[](index);
		return index;
	}

	/**
	 * \brief %Transform a uniformly distributed sample to the stored distribution
	 *
	 * The original sample is value adjusted so that it can be "reused".
	 *
	 * \param[in, out] sampleValue
	 *     An uniformly distributed sample on [0,1]
	 * \return
	 *     The discrete index associated with the sample
	 */
	inline size_t sampleReuse(T &sampleValue) const {
		size_t index = sample(sampleValue);
		sampleValue = (sampleValue - m_cdf[index])
			/ (m_cdf[index + 1] - m_cdf[index]);
		return index;
	}

	/**
	 * \brief %Transform a uniformly distributed sample.
	 *
	 * The original sample is value adjusted so that it can be "reused".
	 *
	 * \param[in,out]
	 *     An uniformly distributed sample on [0,1]
	 * \param[out] pdf
	 *     Probability value of the sample
	 * \return
	 *     The discrete index associated with the sample
	 */
	inline size_t sampleReuse(T &sampleValue, T &pdf) const {
		size_t index = sample(sampleValue, pdf);
		sampleValue = (sampleValue - m_cdf[index])
			/ (m_cdf[index + 1] - m_cdf[index]);
		return index;
	}

//private:
	std::vector<T> m_cdf;
	T m_sum, m_normalization;
	bool m_normalized;
};
typedef DiscretePDF<double> DiscretePDFd;
