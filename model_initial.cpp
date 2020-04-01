#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include "json.hpp"

using json = nlohmann::json;

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
const int NUM_STATES = static_cast<std::underlying_type<PersonState>::type> (
    PersonState::NOCORONA_DEAD) + 1;

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

  int id;
  int category;
  PersonState state = PersonState::SUSCEPTIBLE;
  int days_until_next_state = 0;
  std::vector<int> edges;
  // Needed because person can change state from IMMUNE to NOCORANA_ICU. If
  // patient survives ICU then we need to know does it return to IMMUNE or
  // SUSCEPTIBLE state.
  bool is_immune = false;
};

class Graph;
json simulate(Graph &, json);

class Graph {
  public:
    Graph(json config) {
      // Generate num_persons. Put each person in one of the categories.
      int num_persons = config["num_persons"];
      int last_bound = 0;
      std::vector<int> category_bounds;
      for (auto &x : config["category_ratios"]) {
        last_bound += x.get<int>();
        category_bounds.push_back(x);
      }
      for (int i = 0; i < num_persons; ++i) {
        int x = rand() % category_bounds.back();
        int category = std::upper_bound(
            category_bounds.begin(), category_bounds.end(), x)
              - category_bounds.begin();
        Person person(i, category);
        persons.push_back(person);
      }

      // Add half_degree indirected connections for each person, resulting in
      // 2 * half_degree average degree of a person.
      for (int i = 0; i < num_persons; ++i) {
        auto &person = persons[i];
        for (int j = 0; j < config["half_degree"]; ++j) {
          while (true) {
            int connected_person = rand() % num_persons;
            if (connected_person == person.id) continue;
            if (std::find(person.edges.begin(),
                          person.edges.end(),
                          connected_person) != person.edges.end()) continue;
            person.edges.push_back(connected_person);
            persons[connected_person].edges.push_back(i);
            break;
          }
        }
      }
    }

  private:
    friend json simulate(Graph &, json);
    std::vector<Person> persons;
};

bool bool_with_probability(double p) {
  return rand() < RAND_MAX * p;
}

json simulate(Graph &g, json simulation_config) {
  int num_days = simulation_config["num_days"];
  int num_icus_left = simulation_config["num_icus"];
  auto params_for_categories = simulation_config["initial_params"];
  std::unordered_map<std::string, std::vector<int>> num_per_state_history;
  std::vector<json> events = simulation_config["events"];
  sort(events.begin(), events.end(), [](const json &x, const json &y) {
      return x["day"] < y["day"];
  });
  auto event = events.begin();

  // TODO: Slowest part of the program is probably json lookups in this loop.
  // We can avoid that by extracting data from json to struct outside of the
  // loop.
  for (int i = 0; i < num_days; ++i) {
    std::cerr << "Simulating day " << i << "/" << num_days << "\n";

    while (event != events.end() && event->at("day") == i) {
      auto update_params = event->at("update_params");
      if (update_params.size() != params_for_categories.size()) {
        std::cerr << "illegal number of categories for some event\n";
        exit(1);
      }
      for (int j = 0; j < update_params.size(); ++j) {
        for (auto param : update_params[j].items()) {
          params_for_categories[j][param.key()] = param.value();
        }
      }
      ++event;
    }

    for (auto &x : g.persons) {
      const auto &params = params_for_categories[x.category];

      if (x.state == PersonState::SUSCEPTIBLE) {
        if (bool_with_probability(params["prob_s_to_i"])) {
          // Person is imported case.
          x.state = PersonState::INFECTIOUS;
          x.days_until_next_state = params["days_i_to_c"];
        } else {
          for (int y_id : x.edges) {
            // There are some scenarios that are not clearly defined how they
            // will be resolved, i. e. if y_id < x.id some situation will be
            // treated in one way and if not on the other.
            // For example if person y is infected on this day he can spread
            // infection to person x on the same day if y_id < x.id, otherwise
            // it can't. This problem can be resolved by having state update
            // in multiple passes or by having additional states but probably
            // we don't care much about that.
            const auto &y = g.persons[y_id];
            const auto &y_params = params_for_categories[y.category];
            if (y.state == PersonState::INFECTIOUS &&
                bool_with_probability(y_params["prob_transmission"])) {
              // Person is infected by another infectious person.
              x.state = PersonState::INFECTIOUS;
              x.days_until_next_state = params["days_i_to_c"];
            }
          }
        }

      } else if (x.state == PersonState::INFECTIOUS) {
        if (!--x.days_until_next_state) {
          x.days_until_next_state = params["days_to_im_or_d"];
          // Person became symptomatic so he is either put in isolation or in
          // icu.
          if (bool_with_probability(params["prob_i_to_ic"])) {
            if (!num_icus_left) {
              // Person need icu but there are none left so he dies.
              x.state = PersonState::DEAD;
            } else {
              x.state = PersonState::ICU;
              --num_icus_left;
            }
          } else {
            x.state = PersonState::CONFIRMED;
          }
        }

      } else if (x.state == PersonState::CONFIRMED) {
        if (!--x.days_until_next_state) {
          x.state = PersonState::IMMUNE;
          x.is_immune = true;
        }

      } else if (x.state == PersonState::ICU) {
        if (!--x.days_until_next_state) {
          if (bool_with_probability(params["prob_ic_to_d"])) {
            // Person died in ICU.
            x.state = PersonState::DEAD;
          } else {
            // Person made it, so he became immune.
            x.state = PersonState::IMMUNE;
            x.is_immune = true;
          }
          ++num_icus_left;
        }

      } else if (x.state == PersonState::NOCORONA_ICU) {
        if (!--x.days_until_next_state) {
          if (bool_with_probability(params["prob_nic_to_d"])) {
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

      if (x.state == PersonState::IMMUNE || x.state == PersonState::SUSCEPTIBLE) {
        // Person can require ICU from other illnesses, not only corona.
        if (bool_with_probability(params["prob_to_nic"])) {
          if (num_icus_left) {
            x.state = PersonState::NOCORONA_ICU;
            x.days_until_next_state = params["days_nic"];
            --num_icus_left;
          } else {
            x.state = PersonState::NOCORONA_DEAD;
          }
        }
      }
    }

    // Collecting stats.
    std::unordered_map<std::string, int> num_per_state;
    for (const auto &x : g.persons) {
      ++num_per_state[state_to_name(x.state)];
    }

    // Iterating over all possible states instead of states in num_per_state
    // map since it is not necessary that all possible states are there, but
    // we still want to append zero to history vector.
    for (int j = 0; j < NUM_STATES; ++j) {
      auto state_name = state_to_name(static_cast<PersonState>(j));
      num_per_state_history[state_name].push_back(num_per_state[state_name]);
    }
  }

  return num_per_state_history;
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Expected arguments: config_file seed\n"
              << "  config_file = path to json configuration file\n"
              << "                (see example_config.json for format)\n"
              << "  seed = number passed to srand\n\n";
    exit(1);
  }
  srand(atoi(argv[2]));

  std::ifstream config_stream(argv[1]);
  json config;
  config_stream >> config;
  config_stream.close();

  Graph g(config["graph_generation"]);
  auto stats = simulate(g, config["simulation"]);
  std::cout << stats << "\n";
  return 0;
}
