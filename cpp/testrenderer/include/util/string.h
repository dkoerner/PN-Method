#pragma once
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>



	inline std::string trim_right (const std::string & s, const std::string & t = " \t\r\n");

	inline std::string trim_left (const std::string & s, const std::string & t = " \t\r\n");

	inline std::string trim (const std::string & s, const std::string & t = " \t\r\n");

	std::string replace(std::string &s, const std::string& toReplace, const std::string& replaceWith);
	// split a line into the first word, and rest-of-the-line
	std::string getWord (std::string & s, const std::string delim = " ",const bool trim_spaces = true);

	void splitString( const std::string s, std::vector<std::string> & v, const std::string delim = " ", const bool trim_spaces = true);

	std::string zeroPadNumber(int number, int width = 4);

	std::string toLower( const std::string input );
	bool endsWith(std::string const & value, std::string const & ending);


	// ------------------- FileName / Path Manipulation routines ------------------

	std::string setExtension( const std::string &fileName, const std::string &extension ); // Returns the given filename with another extension.
	std::string getExtension( const std::string &fileName );  // Returns the extension (including period) of the given path

	// ------------------- string - type conversion  ------------------
	template<class T>
    inline std::string toString(const T& t)
	{
		std::ostringstream stream;
		stream << t;
		return stream.str();
	}

	template<class T>
    T fromString(const std::string& s)
	{
		std::istringstream stream (s);
		T t;
		stream >> t;
		return t;
	}
