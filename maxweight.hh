////////////////////////////////////////////////////////////////////////////////
// maxweight.hh
//
// Compute the set of foods that maximizes the weight in foods, within 
// a given maximum calorie amount with the dynamic programming or exhaustive search.
//
///////////////////////////////////////////////////////////////////////////////


#pragma once


#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;

// One food item available for purchase.
class FoodItem
{
	//
	public:
		
		//
		FoodItem
		(
			const std::string& description,
			double calories,
            double weight_ounces
		)
			:
			_description(description),
			_calories(calories),
            _weight_ounces(weight_ounces)
		{
			assert(!description.empty());
			assert(calories > 0);
		}
		
		//
		const std::string& description() const { return _description; }
        double calorie() const { return _calories; }
		double weight() const { return _weight_ounces; }
		
	//
	private:
		
		// Human-readable description of the food, e.g. "spicy chicken breast". Must be non-empty.
		std::string _description;
		
		// Calories; Must be positive
		double _calories;
    
    // Food weight, in ounces; most be non-negative.
		double _weight_ounces;	
};


// Alias for a vector of shared pointers to FoodItem objects.
typedef std::vector<std::shared_ptr<FoodItem>> FoodVector;


// Load all the valid food items from the CSV database
// Food items that are missing fields, or have invalid values, are skipped.
// Returns nullptr on I/O error.
std::unique_ptr<FoodVector> load_food_database(const std::string& path)
{
	std::unique_ptr<FoodVector> failure(nullptr);
	
	std::ifstream f(path);
	if (!f)
	{
		std::cout << "Failed to load food database; Cannot open file: " << path << std::endl;
		return failure;
	}
	
	std::unique_ptr<FoodVector> result(new FoodVector);
	
	size_t line_number = 0;
	for (std::string line; std::getline(f, line); )
	{
		line_number++;
		
		// First line is a header row
		if ( line_number == 1 )
		{
			continue;
		}
		
		std::vector<std::string> fields;
		std::stringstream ss(line);
		
		for (std::string field; std::getline(ss, field, '^'); )
		{
			fields.push_back(field);
		}
		
		if (fields.size() != 3)
		{
			std::cout
				<< "Failed to load food database: Invalid field count at line " << line_number << "; Want 3 but got " << fields.size() << std::endl
				<< "Line: " << line << std::endl
				;
			return failure;
		}
		
		std::string
			descr_field = fields[0],
			calories_field = fields[1],
            weight_ounces_field = fields[2]
			;
		
		auto parse_dbl = [](const std::string& field, double& output)
		{
			std::stringstream ss(field);
			if ( ! ss )
			{
				return false;
			}
			
			ss >> output;
			
			return true;
		};
		
		std::string description(descr_field);
		double calories, weight_ounces;
		if (
			parse_dbl(calories_field, calories)
			&& parse_dbl(weight_ounces_field, weight_ounces)
		)
		{
			result->push_back(
				std::shared_ptr<FoodItem>(
					new FoodItem(
						description,
						calories,
                        weight_ounces
					)
				)
			);
		}
	}

	f.close();
	
	return result;
}


// Convenience function to compute the total weight and calories in 
// a FoodVector.
// Provide the FoodVector as the first argument
// The next two arguments will return the weight and calories back to 
// the caller.
void sum_food_vector
(
	const FoodVector& foods,
	double& total_calories,
    double& total_weight
)
{
	total_calories = total_weight = 0;
	for (auto& food : foods)
	{
		total_calories += food->calorie();
        total_weight += food->weight();
	}
}


// Convenience function to print out each FoodItem in a FoodVector,
// followed by the total weight and calories of it.
void print_food_vector(const FoodVector& foods)
{
	std::cout << "*** food Vector ***" << std::endl;
	
	if ( foods.size() == 0 )
	{
		std::cout << "[empty food list]" << std::endl;
	}
	else
	{
		for (auto& food : foods)
		{
			std::cout
				<< "Ye olde " << food->description()
				<< " ==> "
				<< "; calories = " << food->calorie()
                << "Weight of " << food->weight() << " ounces"
				<< std::endl
				;
		}
		
		double total_calories, total_weight;
		sum_food_vector(foods, total_calories, total_weight);
		std::cout
			<< "> Grand total calories: " << total_calories
			<< std::endl
            << "> Grand total weight: " << total_weight << " ounces" << std::endl
			;
	}
}


// Filter the vector source, i.e. create and return a new FoodVector
// containing the subset of the food items in source that match given
// criteria.
// This is intended to:
//	1) filter out food with zero or negative weight that are irrelevant to // our optimization
//	2) limit the size of inputs to the exhaustive search algorithm since it // will probably be slow.
//
// Each food item that is included must have at minimum min_weight and 
// at most max_weight.
//	(i.e., each included food item's weight must be between min_weight
// and max_weight (inclusive).
//
// In addition, the vector includes only the first total_size food items
// that match these criteria.
std::unique_ptr<FoodVector> filter_food_vector
(
	const FoodVector& source,
	double min_weight,
	double max_weight,
	int total_size
)
{
 // the filtered list and its size
std::unique_ptr<FoodVector> result(new FoodVector);
int size = 0;

// for every FoodItem
for (int index = 0; index < source.size() - 1; index++)
{
	auto foodItemPointer = source[index];
	double weight = foodItemPointer->weight();
	
	if ((weight >= min_weight && weight <= max_weight) && (size < total_size))
	{
		result->push_back(foodItemPointer); // add the food item pointer
		size++;
	}
}

return result;
}


// Compute the optimal set of food items with a exhaustive search algorithm.
// Specifically, among all subsets of food items, return the subset 
// whose weight in ounces fits within the total_weight one can carry and
// whose total calories is greatest.
// To avoid overflow, the size of the food items vector must be less than 64.
/*std::unique_ptr<FoodVector> exhaustive_max_calories
(
	const FoodVector& foods,
	double total_weight
)
{
int n = foods.size();
double bestCalorie = 0.0;
std::unique_ptr<FoodVector> result(new FoodVector); // pointer the the final FoodVector

for (uint64_t bits = 0; bits <= std::pow(2, n) - 1; bits++)
{
	std::unique_ptr<FoodVector> candidate(new FoodVector); // candidate subset
	
	for (int j = 0; j <= n - 1; j++)
	{
		if (((bits >> j) & 1) == 1)
		{
			candidate->push_back(foods[j]);
		}
	}

	// sum up the candidate's weight and calories
	double candidateCalories = 0.0;
	double candidateWeight = 0.0;

	for (auto iter = candidate->begin(); iter != candidate->end(); iter++)
	{
		candidateCalories += (*iter)->calorie();
		candidateWeight += (*iter)->weight();
	}
	
	if (candidateWeight <= total_weight)
		if (result->empty() || candidateCalories > bestCalorie)
		{			
			result->clear();
			
			for (auto iter = candidate->begin(); iter != candidate->end(); iter++)
			{
				result->push_back(*iter);
			}
			
			// sum up the best's calories
			bestCalorie = 0.0;
			for (auto iter = result->begin(); iter != result->end(); iter++)
			{
				bestCalorie += (*iter)->calorie();
			}
		}

}

return result;
}*/

// Compute the optimal set of food items with dynamic programming.
// Specifically, among the food items that fit within a total_calories,
// choose the foods whose weight-per-calorie is largest.
// Repeat until no more food items can be chosen, either because we've 
// run out of food items, or run out of space.
std::unique_ptr<FoodVector> dynamic_max_weight
(
	const FoodVector& foods,
	double total_calories
)
{
std::vector<std::vector<double>> cache;
std::unique_ptr<FoodVector> items(new FoodVector(foods));
std::unique_ptr<FoodVector> best(new FoodVector); // list of food items used to maximize weight
int maxItemIndex = items->size(); // total items along left column. Correlates to n in pseudocode
int maxCalories = int(total_calories); // total calories along top row. Correlates to w in pseudocode
std::vector<double> weights = { 0 }; // weights (values). Leave [0] as 0 because indexing for pseudocode starts at 1. Correlates to V in pseudocode
std::vector<int> calories = { 0 }; // calories (costs). Leave [0] as 0 because indexing for pseudocode starts at 1. Correlates to X in pseudocode

// initialize vectors weights and calories
for (const auto& element : foods)
{
	weights.push_back(element->weight());
	calories.push_back(element->calorie());
}

// set cache dimensions
cache.resize(maxItemIndex + 1);
for (auto& element : cache)
	element.assign(maxCalories + 1, 0); // fill with 0's


/**************** Populate the cache ****************/
for (int i = 1; i <= maxItemIndex; i++) // for items {1, 2, ... maxItemIndex}
	for (int j = 0; j <= maxCalories; j++) // for calories {0, 1, ... maxCalories}
	{
		if (j - calories[i] < 0) // avoid out-of-bounds access
			cache[i][j] = cache[i - 1][j];
		else
			if (cache[i - 1][j - calories[i]] + weights[i] > cache[i - 1][j])
				cache[i][j] = cache[i - 1][j - calories[i]] + weights[i];
			else
				cache[i][j] = cache[i - 1][j];
		
		// display computation of cache[i][j]
		/*if (maxCalories <= 14)
		{
			//cout << "\n\ncalories[" << i << "] = " << calories[i];
			//cout << ", weights[" << i << "] = " << weights[i];
			cout << "\ncache[" << i << "][" << j << "] = max(cache[" << i << "-1][" 
			     << j << "], cache[" << i <<"-1][" << j << "-calories[" << i << "]] + weights[" << i << "]) = " << cache[i][j];
		}*/
	}
// display maximum weight possible within a calorie limit
// cout << "\ncache[maxItemIndex][maxCalories]: " << cache[maxItemIndex][maxCalories];

// rebuild the list of items used
int item = maxItemIndex;
int calorieLimit = maxCalories;
while (item != 0 && calorieLimit != 0)
{
	// display the two options for backtracking
	/*if (maxCalories <= 14)
	{
		cout << "\ncache[item (" << item << ")][calorieLimit (" << calorieLimit << ")]: " << cache[item][calorieLimit];
		cout << "\ncache[item-1][calorieLimit]: " << cache[item-1][calorieLimit];
	}*/
	
	if (cache[item][calorieLimit] == cache[item - 1][calorieLimit]) // if the cell is the same as above
		item--; // go to the row above
	else // if the cell is different from the cell above
	{
		// display the item being included
		/*if (maxCalories <= 14)
			cout << "\nIncluding item: " << item << "\n";*/
		
		// include the item, go X[item] places to the left and up a row
		best->push_back(foods[item - 1]);
		calorieLimit -= calories[item];
		item--;
	}
}
return best;
}

// Compute the optimal set of food items with a exhaustive search algorithm.
// Specifically, among all subsets of food items, return the subset 
// whose weight in ounces fits within the total_weight one can carry and
// whose total calories is greatest.
// To avoid overflow, the size of the food items vector must be less than 64.
std::unique_ptr<FoodVector> exhaustive_max_weight
(
	const FoodVector& foods,
	double total_calorie
)
{
int n = foods.size();
double bestWeight = 0.0;
std::unique_ptr<FoodVector> result(new FoodVector); // pointer the the final FoodVector

for (uint64_t bits = 0; bits <= std::pow(2, n) - 1; bits++)
{
	std::unique_ptr<FoodVector> candidate(new FoodVector); // candidate subset
	
	for (int j = 0; j <= n - 1; j++)
	{
		if (((bits >> j) & 1) == 1)
		{
			candidate->push_back(foods[j]);
		}
	}

	// sum up the candidate's weight and calories
	double candidateCalories = 0.0;
	double candidateWeight = 0.0;

	for (auto iter = candidate->begin(); iter != candidate->end(); iter++)
	{
		candidateCalories += (*iter)->calorie();
		candidateWeight += (*iter)->weight();
	}
	
	if (candidateCalories <= total_calorie)
		if (result->empty() || candidateWeight > bestWeight)
		{			
			result->clear();
			
			for (auto iter = candidate->begin(); iter != candidate->end(); iter++)
			{
				result->push_back(*iter);
			}
			
			// sum up the best's weight
			bestWeight = 0.0;
			for (auto iter = result->begin(); iter != result->end(); iter++)
			{
				bestWeight += (*iter)->weight();
			}
		}

}

return result;
}
