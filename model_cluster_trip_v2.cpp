#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "json.hpp"

using json = nlohmann::json;
using RandomGenerator = std::minstd_rand;

enum class PersonState {
  SUSCEPTIBLE, // s
  INFECTIOUS, // i
  CONFIRMED, // c
  ICU, // ic
  DEAD, // d
  IMMUNE, // im
  NOCORONA_ICU, //nic
  NOCORONA_DEAD, // nd
};
int person_state_to_int(PersonState state) {
  return static_cast<std::underlying_type<PersonState>::type>(state);
}
const int NUM_STATES = person_state_to_int(PersonState::NOCORONA_DEAD) + 1;

std::string state_to_name(PersonState state) {
  switch (state) {
    case PersonState::SUSCEPTIBLE:
      return "susceptible";
    case PersonState::INFECTIOUS:
      return "infectious";
    case PersonState::CONFIRMED:
      return "confirmed";
    case PersonState::ICU:
      return "icu";
    case PersonState::DEAD:
      return "dead";
    case PersonState::IMMUNE:
      return "immune";
    case PersonState::NOCORONA_ICU:
      return "nocorona_icu";
    case PersonState::NOCORONA_DEAD:
      return "nocorona_dead";
  }
  // Can't happen.
  exit(1);
}

struct Person {
  Person(int id, int category) : id(id), category(category) {}
  Person(int id, int category, PersonState state) :
    id(id),
    category(category),
    state(state),
    is_immune(state == PersonState::IMMUNE) {}

  int id;
  int category;
  PersonState state = PersonState::SUSCEPTIBLE;
  int days_until_next_state = 0;
  // Needed because person can change state from IMMUNE to NOCORANA_ICU. If
  // patient survives ICU then we need to know does it return to IMMUNE or
  // SUSCEPTIBLE state.
  bool is_immune = false;
};

struct CategoryParams {
  CategoryParams(json params) :
      prob_goes_on_trip(params["prob_goes_on_trip"]),
      prob_c_trip_candidate(params["prob_c_trip_candidate"]),
      prob_c_neighbour_trip_candidate(
          params["prob_c_neighbour_trip_candidate"]),
      prob_s_to_i(params["prob_s_to_i"]),
      days_i_to_c(params["days_i_to_c"]),
      prob_i_to_ic(params["prob_i_to_ic"]),
      days_c_to_im(params["days_c_to_im"]),
      days_ic_to_im_or_c(params["days_ic_to_im_or_c"]),
      prob_ic_to_d(params["prob_ic_to_d"]),
      prob_to_nic(params["prob_to_nic"]),
      prob_nic_to_d(params["prob_nic_to_d"]),
      days_nic(params["days_nic"]) {}

  double prob_goes_on_trip;
  double prob_c_trip_candidate;
  double prob_c_neighbour_trip_candidate;
  double prob_s_to_i;
  int days_i_to_c;
  double prob_i_to_ic;
  int days_c_to_im;
  int days_ic_to_im_or_c;
  double prob_ic_to_d;
  double prob_to_nic;
  double prob_nic_to_d;
  int days_nic;
};

class Graph;
json simulate(Graph &, json, RandomGenerator &);

class Graph {
  public:
    Graph(json graph_params, RandomGenerator &generator) {
      int person_id = 0;
      for (const auto &subgraph_params: graph_params) {
        // Generate num_persons. Put each person in one of the categories.
        int num_clusters = subgraph_params["num_clusters"];
        int num_people_per_cluster = subgraph_params["num_people_per_cluster"];

        int last_category_bound = 0;
        std::vector<int> category_bounds;
        for (auto &x : subgraph_params["category_ratios"]) {
          last_category_bound += x.get<int>();
          category_bounds.push_back(last_category_bound);
        }

        std::vector<int> people_per_state_ratios =
            subgraph_params.value(
                "people_per_state_ratios",
                std::vector<int>{1, 0, 0, 0, 0, 0, 0, 0});
        if (people_per_state_ratios.size() != NUM_STATES) {
          std::cerr << "Invalid size of people_per_state_ratios\n";
          exit(1);
        }
        int last_state_bound = 0;
        std::vector<int> state_bounds;
        for (int x : people_per_state_ratios) {
          last_state_bound += x;
          state_bounds.push_back(last_state_bound);
        }

        std::uniform_int_distribution<> category_distribution(
            0, category_bounds.back() - 1);
        std::uniform_int_distribution<> state_distribution(
            0, state_bounds.back() - 1);
        for (int i = 0; i < num_clusters; ++i) {
          std::vector<Person> cluster;
          for (int j = 0; j < num_people_per_cluster; ++j) {
            int x = category_distribution(generator);
            int category = std::upper_bound(
                category_bounds.begin(), category_bounds.end(), x)
                  - category_bounds.begin();
            int y = state_distribution(generator);
            auto state = static_cast<PersonState>(std::upper_bound(
                state_bounds.begin(), state_bounds.end(), y)
                  - state_bounds.begin());
            Person person(person_id++, category, state);
            cluster.push_back(person);
          }
          clusters.push_back(cluster);
        }
      }
    }

  private:
    friend json simulate(Graph &, json, RandomGenerator &);
    std::vector<std::vector<Person>> clusters;
};

class BoolWithProbability {
  public:
    BoolWithProbability(RandomGenerator &generator) : generator(generator) {}

    bool operator()(double p) {
      return distribution(generator) < p;
    }

  private:
    RandomGenerator &generator;
    std::uniform_real_distribution<> distribution{0, 1};
};

double dying_probability(double p, double mu, double system_load) {
  return 1 - (1 - p) * exp(-mu * system_load);
}

// Returns whether icu overflow happened.
inline bool before_trip_cluster_update(
    std::vector<Person> &cluster,
    int &num_icus_left,
    const std::vector<CategoryParams> &params_for_categories,
    double prob_transmission,
    double mu,
    double system_load,
    BoolWithProbability &bool_with_probability) {
  // Count number of infected people that can transmit corona virus in this
  // cluster. We do this only once before calculating transitions of people
  // so that all people are in analog position.
  int cnt_infectious_persons = 0;
  for (const auto &x : cluster) {
    if (x.state == PersonState::INFECTIOUS) {
      ++cnt_infectious_persons;
    }
  }
  double p_in_cluster_transmission =
      1 - pow(1 - prob_transmission, cnt_infectious_persons);

  bool icu_overflow = false;
  for (int i = 0; i < cluster.size(); ++i) {
    auto &x = cluster[i];
    const auto &params = params_for_categories[x.category];

    if (x.state == PersonState::SUSCEPTIBLE) {
      if (bool_with_probability(params.prob_s_to_i) ||
          bool_with_probability(p_in_cluster_transmission)) {
        x.state = PersonState::INFECTIOUS;
        x.days_until_next_state = params.days_i_to_c;
      }

    } else if (x.state == PersonState::INFECTIOUS) {
      if (!--x.days_until_next_state) {
        // Person became symptomatic so he is either put in isolation or in
        // icu.
        if (bool_with_probability(params.prob_i_to_ic)) {
          if (!num_icus_left) {
            // Person need icu but there are none left so he dies.
            x.state = PersonState::DEAD;
            icu_overflow = true;
          } else {
            x.state = PersonState::ICU;
            x.days_until_next_state = params.days_ic_to_im_or_c;
            --num_icus_left;
          }
        } else {
          x.state = PersonState::CONFIRMED;
          x.days_until_next_state = params.days_c_to_im;
        }
      }

    } else if (x.state == PersonState::CONFIRMED) {
      if (!--x.days_until_next_state) {
        x.state = PersonState::IMMUNE;
        x.is_immune = true;
      }

    } else if (x.state == PersonState::ICU) {
      if (!--x.days_until_next_state) {
        double p_ic_to_d = dying_probability(
            params.prob_ic_to_d, mu, system_load);
        if (bool_with_probability(p_ic_to_d)) {
          // Person died in ICU.
          x.state = PersonState::DEAD;
        } else {
          // Person made it through ICU, so he is now mild case.
          x.state = PersonState::CONFIRMED;
          x.days_until_next_state = params.days_c_to_im;
        }
        ++num_icus_left;
      }

    } else if (x.state == PersonState::NOCORONA_ICU) {
      if (!--x.days_until_next_state) {
        double p_nic_to_d = dying_probability(
            params.prob_nic_to_d, mu, system_load);
        if (bool_with_probability(p_nic_to_d)) {
          // Person died in ICU of illness that is not corona.
          x.state = PersonState::NOCORONA_DEAD;
        } else if (x.is_immune) {
          x.state = PersonState::IMMUNE;
        } else {
          x.state = PersonState::SUSCEPTIBLE;
        }
        ++num_icus_left;
      }
    }

    if (x.state == PersonState::IMMUNE ||
        x.state == PersonState::SUSCEPTIBLE ||
        x.state == PersonState::CONFIRMED ||
        x.state == PersonState::INFECTIOUS) {
      // Person can require ICU from other illnesses, not only corona.
      if (bool_with_probability(params.prob_to_nic)) {
        if (num_icus_left) {
          if (x.state == PersonState::CONFIRMED ||
              x.state == PersonState::INFECTIOUS) {
            // If person who had corona went to an icu for unrelated reasons we
            // presume that he will get over that corona infection during his
            // time in icu. That does not neccessarily follow from given days
            // for NOCORONA_ICU state and different stages of corona disease
            // from config but simplifies implementation and has almost no
            // impact on simulation since number of NOCORONA_ICU people is
            // expected to be small.
            x.is_immune = true;
          }
          x.state = PersonState::NOCORONA_ICU;
          x.days_until_next_state = params.days_nic;
          --num_icus_left;
        } else {
          x.state = PersonState::NOCORONA_DEAD;
          icu_overflow = true;
        }
      }
    }
  }

  return icu_overflow;
}

json simulate(Graph &g, json simulation_config, RandomGenerator &generator) {
  BoolWithProbability bool_with_probability(generator);
  // Extracting configuration parameters.
  int num_days = simulation_config["stopping_conditions"]["num_days"];
  bool on_icu_overflow =
    simulation_config["stopping_conditions"]["on_icu_overflow"];
  bool on_pandemic_end =
    simulation_config["stopping_conditions"].value("on_pandemic_end", false);
  int num_icus_left = simulation_config["num_icus"];
  double mu = simulation_config["mu"];
  double prob_transmission = simulation_config["prob_transmission"];
  double k_trip = simulation_config["k_trip"];
  bool isolate_cluster_on_known_case =
      simulation_config["isolate_cluster_on_known_case"];

  std::vector<CategoryParams> params_for_categories;
  auto all_params = simulation_config["initial_params"];
  for (const auto &params : all_params) {
    params_for_categories.emplace_back(params);
  }

  std::vector<json> events = simulation_config["events"];
  sort(events.begin(), events.end(), [](const json &x, const json &y) {
      return x["day"] < y["day"];
  });
  auto event = events.begin();

  // Declaring stats variables.
  std::unordered_map<std::string, std::vector<int>> num_per_state_history;
  std::string stopping_condition = "num_days";
  int num_days_icu_overflow = 0;
  int first_day_icu_overflow = -1;
  int last_day_icu_overflow = -1;

  for (int day = 0; day < num_days; ++day) {
    std::cerr << "Simulating day " << day << "/" << num_days << "\n";
    bool this_day_icu_overflow = false;

    while (event != events.end() && event->at("day") == day) {
      auto update_params = event->at("update_params");
      for (auto param: update_params.items()) {
        if (!all_params[0].count(param.key())) {
          std::cerr << "Invalid key `" << param.key() << "` in an event\n";
          exit(1);
        }
        if (param.value().size() != all_params.size()) {
          std::cerr << "Invalid number of categories for key `" << param.key()
            << "` in an event\n";
          exit(1);
        }
        for (int i = 0; i < all_params.size(); ++i) {
          all_params[i][param.key()] = param.value()[i];
        }
      }
      params_for_categories.clear();
      for (const auto &params : all_params) {
        params_for_categories.emplace_back(params);
      }
      ++event;
    }

    // Calculate system_load factor.
    int cnt_alive_people = 0;
    int cnt_burden = 0;
    for (const auto &cluster : g.clusters) {
      for (const auto &x : cluster) {
        if (x.state != PersonState::DEAD &&
            x.state != PersonState::NOCORONA_DEAD) {
          ++cnt_alive_people;
        }
        if (x.state == PersonState::ICU ||
            x.state == PersonState::CONFIRMED) {
          ++cnt_burden;
        }
      }
    }
    double system_load = (double)cnt_burden / cnt_alive_people;

    // Before "trip" updates.
    for (auto &cluster : g.clusters) {
      bool cluster_icu_overflow = before_trip_cluster_update(
          cluster,
          num_icus_left,
          params_for_categories,
          prob_transmission,
          mu,
          system_load,
          bool_with_probability);
      if (cluster_icu_overflow && !this_day_icu_overflow) {
        this_day_icu_overflow = true;
        ++num_days_icu_overflow;
        last_day_icu_overflow = day;
        if (first_day_icu_overflow == -1) {
          first_day_icu_overflow = day;
        }
      }
    }

    // Pick people who go to the trip.
    std::vector<Person*> persons_on_trip;
    int num_people_on_trip_with_cluster_corona = 0;
    int num_able_people_with_cluster_corona = 0;
    for (auto &cluster : g.clusters) {
      // Check if there is someone with known corona disease in cluster.
      bool has_known_corona = false;
      for (auto &x : cluster) {
        if (x.state == PersonState::ICU ||
            x.state == PersonState::CONFIRMED) {
          // It can happen that person gets to NOCORONA_ICU and already had
          // corona. In that case this flag would stay false, however I don't
          // think it is a problem since this should happen quite rarely.
          has_known_corona = true;
        }
      }

      for (auto &x : cluster) {
        const auto &params = params_for_categories[x.category];
        if (x.state == PersonState::SUSCEPTIBLE ||
            x.state == PersonState::INFECTIOUS ||
            x.state == PersonState::IMMUNE) {
          if (has_known_corona) {
            ++num_able_people_with_cluster_corona;
          }
          if (!has_known_corona || !isolate_cluster_on_known_case ||
              bool_with_probability(params.prob_c_neighbour_trip_candidate)) {
            if (bool_with_probability(params.prob_goes_on_trip)) {
              persons_on_trip.push_back(&x);
              if (has_known_corona) {
                ++num_people_on_trip_with_cluster_corona;
              }
            }
          }
        }
        if (x.state == PersonState::CONFIRMED) {
          // Person knows that it has corona but it can disobey order for
          // staying home and becomes trip candidate.
          if (bool_with_probability(
                params.prob_c_trip_candidate * params.prob_goes_on_trip)) {
            persons_on_trip.push_back(&x);
            ++num_people_on_trip_with_cluster_corona;
          }
        }
      }
    }

    // Count number of infectious and confirmed people on a trip.
    int cnt_contagious = 0;
    for (Person *x: persons_on_trip) {
      if (x->state == PersonState::INFECTIOUS ||
          x->state == PersonState::CONFIRMED) {
        ++cnt_contagious;
      }
    }

    // Spread infection during the trip.
    double contagious_ratio = (double)cnt_contagious / persons_on_trip.size();
    double p_transmission = std::min(
        prob_transmission * k_trip * contagious_ratio, 1.0);
    for (Person *x: persons_on_trip) {
      if (x->state == PersonState::SUSCEPTIBLE &&
          bool_with_probability(p_transmission)) {
        const auto &params = params_for_categories[x->category];
        x->state = PersonState::INFECTIOUS;
        x->days_until_next_state = params.days_i_to_c;
      }
    }

    // Collecting stats.
    std::vector<int> num_per_state(NUM_STATES);
    for (const auto &cluster : g.clusters) {
      for (const auto &x : cluster) {
        ++num_per_state[person_state_to_int(x.state)];
      }
    }

    // Iterating over all possible states instead of states in num_per_state
    // map since it is not necessary that all possible states are there, but
    // we still want to append zero to history vector.
    for (int j = 0; j < NUM_STATES; ++j) {
      auto state_name = state_to_name(static_cast<PersonState>(j));
      num_per_state_history[state_name].push_back(num_per_state[j]);
    }
    num_per_state_history["num_people_on_trip"].push_back(persons_on_trip.size());
    num_per_state_history["num_people_on_trip_with_cluster_corona"].push_back(
        num_people_on_trip_with_cluster_corona);
    num_per_state_history["num_able_people_with_cluster_corona"].push_back(
        num_able_people_with_cluster_corona);

    if (this_day_icu_overflow && on_icu_overflow) {
      stopping_condition = "icu_overflow";
      break;
    }
    if (event == events.end() && on_pandemic_end &&
        num_per_state[person_state_to_int(PersonState::INFECTIOUS)] == 0 &&
        num_per_state[person_state_to_int(PersonState::CONFIRMED)] == 0 &&
        num_per_state[person_state_to_int(PersonState::ICU)] == 0) {
      stopping_condition = "pandemic_end";
      break;
    }
  }

  return {
    {"stopping_condition", stopping_condition},
    {"num_days_icu_overflow", num_days_icu_overflow},
    {"first_day_icu_overflow", first_day_icu_overflow},
    {"last_day_icu_overflow", last_day_icu_overflow},
    {"stats", num_per_state_history},
  };
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Expected arguments: config_file seed\n"
              << "  config_file = path to json configuration file\n"
              << "                (see example_config.json for format)\n"
              << "  seed = number passed to generator constructor\n\n";
    exit(1);
  }

  std::ifstream config_stream(argv[1]);
  json config;
  config_stream >> config;
  config_stream.close();

  RandomGenerator generator(atoi(argv[2]));
  Graph g(config["graph_generation"], generator);
  auto data = simulate(g, config["simulation"], generator);
  data["config"] = config;
  std::cout << data << "\n";
  return 0;
}
