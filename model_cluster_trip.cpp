#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
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
  // Needed because person can change state from IMMUNE to NOCORANA_ICU. If
  // patient survives ICU then we need to know does it return to IMMUNE or
  // SUSCEPTIBLE state.
  bool is_immune = false;
};

struct CategoryParams {
  CategoryParams(json params) :
      prob_c_trip_candidate(params["prob_c_trip_candidate"]),
      prob_s_to_i(params["prob_s_to_i"]),
      days_i_to_c(params["days_i_to_c"]),
      prob_i_to_ic(params["prob_i_to_ic"]),
      days_c_to_im(params["days_c_to_im"]),
      days_ic_to_im_or_d(params["days_ic_to_im_or_d"]),
      prob_ic_to_d(params["prob_ic_to_d"]),
      prob_to_nic(params["prob_to_nic"]),
      prob_nic_to_d(params["prob_nic_to_d"]),
      days_nic(params["days_nic"]) {}

  double prob_c_trip_candidate;
  double prob_s_to_i;
  int days_i_to_c;
  double prob_i_to_ic;
  int days_c_to_im;
  int days_ic_to_im_or_d;
  double prob_ic_to_d;
  double prob_to_nic;
  double prob_nic_to_d;
  int days_nic;
};

class Graph;
json simulate(Graph &, json, std::mt19937 &);

class Graph {
  public:
    Graph(json config, std::mt19937 &generator) {
      // Generate num_persons. Put each person in one of the categories.
      int num_clusters = config["num_clusters"];
      int num_people_per_cluster = config["num_people_per_cluster"];
      int last_bound = 0;
      std::vector<int> category_bounds;
      for (auto &x : config["category_ratios"]) {
        last_bound += x.get<int>();
        category_bounds.push_back(x);
      }

      std::uniform_int_distribution<> distribution(
          0, category_bounds.back() - 1);
      for (int i = 0; i < num_clusters; ++i) {
        std::vector<Person> cluster;
        for (int j = 0; j < num_people_per_cluster; ++j) {
          int x = distribution(generator);
          int category = std::upper_bound(
              category_bounds.begin(), category_bounds.end(), x)
                - category_bounds.begin();
          int id = i * num_people_per_cluster + j;
          Person person(id, category);
          cluster.push_back(person);
        }
        clusters.push_back(cluster);
      }
    }

  private:
    friend json simulate(Graph &, json, std::mt19937 &);
    std::vector<std::vector<Person>> clusters;
};

class BoolWithProbability {
  public:
    BoolWithProbability(std::mt19937 &generator) : generator(generator) {}

    bool operator()(double p) {
      return distribution(generator) < p;
    }

  private:
    std::mt19937 &generator;
    std::uniform_real_distribution<> distribution{0, 1};
};

double dying_probability(double p, double mu, double system_load) {
  return 1 - (1 - p) * exp(-mu * system_load);
}

void before_trip_cluster_update(
    std::vector<Person> &cluster,
    int &num_icus_left,
    const std::vector<CategoryParams> &params_for_categories,
    double prob_transmission,
    double mu,
    double system_load,
    BoolWithProbability &bool_with_probability) {
  for (int i = 0; i < cluster.size(); ++i) {
    auto &x = cluster[i];
    const auto &params = params_for_categories[x.category];

    if (x.state == PersonState::SUSCEPTIBLE) {
      if (bool_with_probability(params.prob_s_to_i)) {
        // Person is imported case.
        x.state = PersonState::INFECTIOUS;
        x.days_until_next_state = params.days_i_to_c;
      } else {
        for (int j = 0; j < cluster.size(); ++j) {
          if (i == j) continue;
          // There are some scenarios that are not clearly defined how they
          // will be resolved, i. e. if y.id < x.id some situation will be
          // treated in one way and if not on the other.
          // For example if person y is infected on this day he can spread
          // infection to person x on the same day if y.id < x.id, otherwise
          // it can't. This problem can be resolved by having state update
          // in multiple passes or by having additional states but probably
          // we don't care much about that.
          const auto &y = cluster[j];
          if (y.state == PersonState::INFECTIOUS &&
              bool_with_probability(prob_transmission)) {
            // Person is infected by another infectious person.
            x.state = PersonState::INFECTIOUS;
            x.days_until_next_state = params.days_i_to_c;
          }
        }
      }

    } else if (x.state == PersonState::INFECTIOUS) {
      if (!--x.days_until_next_state) {
        // Person became symptomatic so he is either put in isolation or in
        // icu.
        if (bool_with_probability(params.prob_i_to_ic)) {
          if (!num_icus_left) {
            // Person need icu but there are none left so he dies.
            x.state = PersonState::DEAD;
          } else {
            x.state = PersonState::ICU;
            x.days_until_next_state = params.days_ic_to_im_or_d;
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
          // Person made it, so he became immune.
          x.state = PersonState::IMMUNE;
          x.is_immune = true;
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
        }
      }
    }
  }
}

json simulate(Graph &g, json simulation_config, std::mt19937 &generator) {
  BoolWithProbability bool_with_probability(generator);
  int num_days = simulation_config["num_days"];
  int num_icus_left = simulation_config["num_icus"];
  double mu = simulation_config["mu"];
  double p_goes_on_trip = simulation_config["prob_goes_on_trip"];
  double prob_transmission = simulation_config["prob_transmission"];
  double k_trip = simulation_config["k_trip"];
  bool isolate_cluster_on_known_case =
      simulation_config["isolate_cluster_on_known_case"];

  std::vector<CategoryParams> params_for_categories;
  for (const auto &params : simulation_config["initial_params"]) {
    params_for_categories.emplace_back(params);
  }

  std::unordered_map<std::string, std::vector<int>> num_per_state_history;
  std::vector<json> events = simulation_config["events"];
  sort(events.begin(), events.end(), [](const json &x, const json &y) {
      return x["day"] < y["day"];
  });
  auto event = events.begin();

  for (int day = 0; day < num_days; ++day) {
    std::cerr << "Simulating day " << day << "/" << num_days << "\n";

    while (event != events.end() && event->at("day") == day) {
      if (event->contains("prob_goes_on_trip")) {
        p_goes_on_trip = event->at("prob_goes_on_trip");
      } else if (event->contains("prob_s_to_i")) {
        auto new_prob_s_to_i = event->at("prob_s_to_i");
        if (new_prob_s_to_i.size() != params_for_categories.size()) {
          std::cerr << "wrong number of categories in event\n";
          exit(1);
        }
        for (int i = 0; i < new_prob_s_to_i.size(); ++i) {
          params_for_categories[i].prob_s_to_i = new_prob_s_to_i[i];
        }
      } else {
        std::cerr << "unsupported param in event\n";
        exit(1);
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
      before_trip_cluster_update(
          cluster,
          num_icus_left,
          params_for_categories,
          prob_transmission,
          mu,
          system_load,
          bool_with_probability);
    }

    // Pick people who go to the trip.
    std::vector<Person*> persons_on_trip;
    for (auto &cluster : g.clusters) {
      // Check if there is someone with known corona disease in cluster.
      bool has_known_corona = false;
      for (auto &x : cluster) {
        if (x.state == PersonState::ICU ||
            x.state == PersonState::CONFIRMED) {
          has_known_corona = true;
        }
      }

      std::vector<Person*> candidates;
      for (auto &x : cluster) {
        const auto &params = params_for_categories[x.category];
        if (x.state == PersonState::SUSCEPTIBLE ||
            x.state == PersonState::INFECTIOUS ||
            x.state == PersonState::IMMUNE) {
          if (!has_known_corona || !isolate_cluster_on_known_case ||
                bool_with_probability(params.prob_c_trip_candidate)) {
            candidates.push_back(&x);
          }
        }
        if (x.state == PersonState::CONFIRMED) {
          // Person knows that it has corona but it can disobey order for
          // staying home and becomes trip candidate.
          if (bool_with_probability(params.prob_c_trip_candidate)) {
            candidates.push_back(&x);
          }
        }
      }
      if (!candidates.empty() && bool_with_probability(p_goes_on_trip)) {
        std::uniform_int_distribution<> distribution(
            0, candidates.size() - 1);
        int candidate_idx = distribution(generator);
        persons_on_trip.push_back(candidates[candidate_idx]);
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
        ++num_per_state[
          static_cast<std::underlying_type<PersonState>::type>(x.state)];
      }
    }

    // Iterating over all possible states instead of states in num_per_state
    // map since it is not necessary that all possible states are there, but
    // we still want to append zero to history vector.
    for (int j = 0; j < NUM_STATES; ++j) {
      auto state_name = state_to_name(static_cast<PersonState>(j));
      num_per_state_history[state_name].push_back(num_per_state[j]);
    }
  }

  return num_per_state_history;
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

  std::mt19937 generator(atoi(argv[2]));
  Graph g(config["graph_generation"], generator);
  auto stats = simulate(g, config["simulation"], generator);
  std::cout << stats << "\n";
  return 0;
}
