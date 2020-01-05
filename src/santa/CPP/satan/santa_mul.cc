#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <string>
#include <queue>
#include <map>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <set>
#include <stdint.h>
#include <sys/time.h>
#include <thread>
#include <mutex>
using namespace std;

class CRandomMother {
public:
   void RandomInit(int seed);
   double Random();
   uint32_t BRandom();
   CRandomMother(int seed) {
      RandomInit(seed);
   }
protected:
   uint32_t x[5];
};

uint32_t CRandomMother::BRandom() {
  uint64_t sum;
  sum = (uint64_t)2111111111UL * (uint64_t)x[3] +
     (uint64_t)1492 * (uint64_t)(x[2]) +
     (uint64_t)1776 * (uint64_t)(x[1]) +
     (uint64_t)5115 * (uint64_t)(x[0]) +
     (uint64_t)x[4];
  x[3] = x[2];  x[2] = x[1];  x[1] = x[0];
  x[4] = (uint32_t)(sum >> 32);                  // Carry
  x[0] = (uint32_t)sum;                          // Low 32 bits of sum
  return x[0];
}

double CRandomMother::Random() {
  return BRandom() / (65536.0 * 65536.0);
}

void CRandomMother::RandomInit (int seed) {
  int i;
  uint32_t s = seed;
  for (i = 0; i < 5; i++) {
    s = s * 29943829 - 1;
    x[i] = s;
  }
  for (i=0; i<19; i++) BRandom();
}

struct _gift_info {
  double lat, lng, weight;
  double x, y, z;
};

struct _trip_info {
  double weight_sum;
  double cur_penalty;
  double min_neighbor;
  vector<int> seq;
  vector<double> acc_weight;
  vector<double> acc_distance;
  vector<double> acc_penalty;
  vector<double> step_distance;
  mutex trip_mutex;
};

const int gift_cnt = 100000;
const int data_cnt = 1000;
const int data_used_bit = 6;

int db_by_dist[gift_cnt][1 << data_used_bit];

_gift_info gifts[gift_cnt+1];
int cur_trip[gift_cnt];

/*double cur_penalty[gift_cnt];
vector<int> trips[100000];
vector<double> acc_weight[100000];
vector<double> acc_distance[100000];
vector<double> acc_penalty[100000];
vector<double> step_distance[100000];
double weight_sum[gift_cnt];*/
_trip_info trips[gift_cnt+1];

const int max_thread_cnt = 100;
double new_change;
long long thread_iter[max_thread_cnt];
double gain[max_thread_cnt];

const double PI = 3.1415926535897932384626;

double getTime()
{
  timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + tv.tv_usec * 0.000001;
}

// temporary verifying function
bool verify(double a, double b) {
  if (fabs(a-b) > 1e-9) {
    printf("! %.14f %.14f %.14f\n", a, b, a-b);
    while (1);
  }
}

double dist(int g1, int g2) {
  if (g1 == gift_cnt) return PI/2 - gifts[g2].lat;
  if (g2 == gift_cnt) return PI/2 - gifts[g1].lat;

  double dx = gifts[g1].x-gifts[g2].x;
  double dy = gifts[g1].y-gifts[g2].y;
  double dz = gifts[g1].z-gifts[g2].z;
  double d = dx*dx+dy*dy+dz*dz;
  double r = sqrt(d);
  if (r < 0.02) return r * (2 + (2./6)*d);
  return 2 * asin(r);
//  double y = gifts[g1].sin_lat*gifts[g2].sin_lat + gifts[g1].cos_lat*gifts[g2].cos_lat*(gifts[g1].sin_lng*gifts[g2].sin_lng+gifts[g1].cos_lng*gifts[g2].cos_lng);
//  double y = sin_lat[g1]*sin_lat[g2] + cos_lat[g1]*cos_lat[g2]*(sin_lng[g1]*sin_lng[g2]+cos_lng[g1]*cos_lng[g2]);
//  if (y > 0.996) {
//    double m0 = 6-sqrt(12+24*y);
//    double m1 = sqrt(m0);
//    return m1 / 720 * (720 - m0*m0);
//  }

//  return acos(y);
//  return acos(sin_lat[g1]*sin_lat[g2] + cos_lat[g1]*cos_lat[g2]*cos(lng[g2]-lng[g1]));
//  double right = cos(lat[g2]-lat[g1]) - cos_lat[g1]*cos_lat[g2]*(1-cos(lng[g2]-lng[g1]));
//  return acos(right);
}

double calc_penalty(const vector<int>& trip, double weight_sum) {
  if (!trip.size()) return 0.0;
  double res = 0.0;
  double w = weight_sum + 10.0;
  for (int i=0; i<trip.size()-1; i++) {
    res += dist(trip[i], trip[i+1]) * w;
    w -= gifts[trip[i+1]].weight;
  }
  return res;
}

int update_trip(int i, int update_from) {
  int sz = trips[i].seq.size() - 2;
  trips[i].acc_weight.resize(sz+1);
  trips[i].acc_distance.resize(sz+2);
  trips[i].acc_penalty.resize(sz+2);
  trips[i].step_distance.resize(sz+1);
  trips[i].acc_weight[0] = trips[i].weight_sum + 10;
  int last_trip;
  if (update_from == 0) last_trip = gift_cnt;
  else last_trip = trips[i].seq[update_from];
  for (int j=update_from; j<sz; j++) {
    cur_trip[trips[i].seq[j+1]] = (i<<11) + j;
    double d = dist(last_trip, trips[i].seq[j+1]);
    trips[i].step_distance[j] = d;
    trips[i].acc_weight[j+1] = trips[i].acc_weight[j] - gifts[trips[i].seq[j+1]].weight;
    trips[i].acc_distance[j+1] = trips[i].acc_distance[j] + d;
    trips[i].acc_penalty[j+1] = trips[i].acc_penalty[j] + trips[i].acc_weight[j] * d;
    last_trip = trips[i].seq[j+1];
  }
  double d = dist(last_trip, gift_cnt);
  trips[i].step_distance[sz] = d;
  trips[i].acc_distance[sz+1] = trips[i].acc_distance[sz] + d;
  trips[i].cur_penalty = trips[i].acc_penalty[sz] + 10 * d;
  trips[i].acc_penalty[sz+1] = trips[i].cur_penalty;
}

int update_trip_without_distance(int i, int update_from, int update_distance_from) {
  int sz = trips[i].seq.size() - 2;
  trips[i].acc_weight.resize(sz+1);
  trips[i].acc_distance.resize(sz+2);
  trips[i].acc_penalty.resize(sz+2);
  trips[i].acc_weight[0] = trips[i].weight_sum + 10;
  for (int j=update_from; j<sz; j++) {
    double d = trips[i].step_distance[j];
    trips[i].acc_weight[j+1] = trips[i].acc_weight[j] - gifts[trips[i].seq[j+1]].weight;
    trips[i].acc_penalty[j+1] = trips[i].acc_penalty[j] + trips[i].acc_weight[j] * d;
  }
  for (int j=update_distance_from; j<sz; j++) {
    cur_trip[trips[i].seq[j+1]] = (i<<11) + j;
    double d = trips[i].step_distance[j];
    trips[i].acc_distance[j+1] = trips[i].acc_distance[j] + d;
  }
  double d = trips[i].step_distance[sz];
  trips[i].acc_distance[sz+1] = trips[i].acc_distance[sz] + d;
  trips[i].cur_penalty = trips[i].acc_penalty[sz] + 10 * d;
  trips[i].acc_penalty[sz+1] = trips[i].cur_penalty;
}

// at gift #idx,
// weight from prev = acc_weight[idx]
// distance to prev = acc_distance[idx+1] - acc_distance[idx]
// penalty to prev = acc_penalty[idx+1] - acc_penalty[idx]

// gift #new_gift is added just before trips[trip_num][idx]
double delta_inserted_gift(int trip_num, int idx, int new_gift, double ww) {
  int gift_prev = trips[trip_num].seq[idx];
  int gift_next = trips[trip_num].seq[idx + 1];

  double res = 0.0;
  res += trips[trip_num].acc_distance[idx] * ww;
  res += dist(gift_prev, new_gift) * (trips[trip_num].acc_weight[idx] + ww);
  res += (dist(new_gift, gift_next) - trips[trip_num].step_distance[idx]) * trips[trip_num].acc_weight[idx];
  return res;
}

// gift trips[trip_num][idx] is removed
double delta_removed_gift(int trip_num, int idx, double ww) {
  int gift_prev = trips[trip_num].seq[idx];
  int gift_next = trips[trip_num].seq[idx + 2];

  double res = 0.0;

  res -= trips[trip_num].acc_distance[idx] * ww;

  res -= trips[trip_num].acc_penalty[idx+2] - trips[trip_num].acc_penalty[idx];
  res += dist(gift_prev, gift_next) * trips[trip_num].acc_weight[idx+1];
  return res;
}

// gift #new_gift is changed from trips[trip_num][idx]
double delta_changed_gift(int trip_num, int idx, int new_gift) {
  int gift_prev = trips[trip_num].seq[idx];
  int gift_next = trips[trip_num].seq[idx + 2];

  double res = 0.0;
  res -= trips[trip_num].acc_distance[idx] * (trips[trip_num].acc_weight[idx] - trips[trip_num].acc_weight[idx+1] - gifts[new_gift].weight);
  res -= trips[trip_num].acc_penalty[idx+2] - trips[trip_num].acc_penalty[idx];
  res += dist(gift_prev, new_gift) * (trips[trip_num].acc_weight[idx+1] + gifts[new_gift].weight);
  res += dist(new_gift, gift_next) * trips[trip_num].acc_weight[idx+1];
  return res;
}


double delta_interchanged_gift(int trip_num, int insert_pos, int new_gift, int remove_pos) {
  double res = 0.0;
  if (remove_pos == insert_pos || remove_pos == insert_pos - 1)
    return delta_changed_gift(trip_num, remove_pos, new_gift);

  if (remove_pos > insert_pos) {
    res = delta_removed_gift(trip_num, remove_pos, gifts[trips[trip_num].seq[remove_pos+1]].weight);
    int gift_prev = trips[trip_num].seq[insert_pos];
    int gift_next = trips[trip_num].seq[insert_pos + 1];

    res += trips[trip_num].acc_distance[insert_pos] * gifts[new_gift].weight;
    res += dist(gift_prev, new_gift) * (trips[trip_num].acc_weight[insert_pos] - gifts[trips[trip_num].seq[remove_pos+1]].weight + gifts[new_gift].weight);
    res += (dist(new_gift, gift_next) - trips[trip_num].step_distance[insert_pos]) * (trips[trip_num].acc_weight[insert_pos] - gifts[trips[trip_num].seq[remove_pos+1]].weight);
  } else {
    res = delta_inserted_gift(trip_num, insert_pos, new_gift, gifts[new_gift].weight);
    int gift_prev = trips[trip_num].seq[remove_pos];
    int gift_next = trips[trip_num].seq[remove_pos + 2];

    res -= trips[trip_num].acc_distance[remove_pos] * gifts[trips[trip_num].seq[remove_pos+1]].weight;

    res -= trips[trip_num].acc_penalty[remove_pos+2] - trips[trip_num].acc_penalty[remove_pos];
    res -= (trips[trip_num].acc_distance[remove_pos+2] - trips[trip_num].acc_distance[remove_pos]) * gifts[new_gift].weight;
    res += dist(gift_prev, gift_next) * (trips[trip_num].acc_weight[remove_pos+1] + gifts[new_gift].weight);
  }
  return res;
}

// trip t1{s1..f1} is switched to t2{s2..f2}
double delta_changed_trip(int t1, int s1, int f1, int t2, int s2, int f2) {
//  printf("%d %d %d %d %d %d\n", t1, s1, f1, t2, s2, f2);
  if (s1 > f1 && s2 > f2) return 0.0;

  int prev1, next1;
  prev1 = trips[t1].seq[s1];
  next1 = trips[t1].seq[f1+2];

  double weight_s1_f1 = -trips[t1].acc_weight[f1+1]+trips[t1].acc_weight[s1];
  double weight_s2_f2 = -trips[t2].acc_weight[f2+1]+trips[t2].acc_weight[s2];
  double res = 0.0;
  double diff = weight_s1_f1 - weight_s2_f2;
  res -= trips[t1].acc_distance[s1] * diff;

  res -= trips[t1].acc_penalty[f1+1] - trips[t1].acc_penalty[s1];
  if (s2 <= f2) {
    res += dist(prev1, trips[t2].seq[s2+1]) * (trips[t1].acc_weight[s1] - diff);
    res += trips[t2].acc_penalty[f2+1] - trips[t2].acc_penalty[s2+1];
    res -= (trips[t2].acc_weight[f2+1] - trips[t1].acc_weight[f1+1]) * (trips[t2].acc_distance[f2+1] - trips[t2].acc_distance[s2+1]);
  }

  if (s2 <= f2) {
    res += (dist(trips[t2].seq[f2+1], next1) - trips[t1].step_distance[f1+1]) * trips[t1].acc_weight[f1+1];
  } else {
    res += (dist(prev1, next1) - trips[t1].step_distance[f1+1]) * trips[t1].acc_weight[f1+1];
  }
  return res;
}

// i1 is moved just before i2
double delta_moved_gift_same_trip(int trip_num, int i1, int i2) {
  int p1, n1, p2, n2;
  p1 = trips[trip_num].seq[i1];
  n1 = trips[trip_num].seq[i1+2];
  p2 = trips[trip_num].seq[i2];
  n2 = trips[trip_num].seq[i2+1];

  double res = 0.0;
  res -= trips[trip_num].acc_penalty[i1 + 2] - trips[trip_num].acc_penalty[i1];
  res -= trips[trip_num].acc_penalty[i2 + 1] - trips[trip_num].acc_penalty[i2];
  double i1_weight = trips[trip_num].acc_weight[i1] - trips[trip_num].acc_weight[i1+1];
  if (i1 < i2) {
    res += dist(p1, n1) * trips[trip_num].acc_weight[i1];
    res += (trips[trip_num].acc_distance[i2] - trips[trip_num].acc_distance[i1 + 2]) * i1_weight;
    res += dist(p2, trips[trip_num].seq[i1+1]) * (trips[trip_num].acc_weight[i2] + i1_weight);
    res += dist(trips[trip_num].seq[i1+1], n2) * trips[trip_num].acc_weight[i2];
  } else {
    res += dist(p1, n1) * trips[trip_num].acc_weight[i1 + 1];
    res -= (trips[trip_num].acc_distance[i1] - trips[trip_num].acc_distance[i2 + 1]) * i1_weight;
    res += dist(p2, trips[trip_num].seq[i1+1]) * trips[trip_num].acc_weight[i2];
    res += dist(trips[trip_num].seq[i1+1], n2) * (trips[trip_num].acc_weight[i2] - i1_weight);
  }
  return res;
}

bool check_delta_swapped_gift(int trip_num, int i1, int i2, double max_margin, double* real_margin) {
  if (i1 > i2) swap(i1, i2);
  if (i1 == i2) return 0.0;

  int gift_prev = trips[trip_num].seq[i1];
  int gift_next = trips[trip_num].seq[i2 + 2];

  double res = 0.0;
  double diff = gifts[trips[trip_num].seq[i1+1]].weight - gifts[trips[trip_num].seq[i2+1]].weight;
  if (i2 == i1 + 1) {
    res += trips[trip_num].step_distance[i2] * diff;
    res -= trips[trip_num].step_distance[i1] * trips[trip_num].acc_weight[i1];
    res -= trips[trip_num].step_distance[i2+1] * trips[trip_num].acc_weight[i2+1];
    if (res > max_margin) return false;
    res += dist(gift_prev, trips[trip_num].seq[i2+1]) * trips[trip_num].acc_weight[i1];
    if (res > max_margin) return false;
    res += dist(trips[trip_num].seq[i1+1], gift_next) * trips[trip_num].acc_weight[i2+1];
  } else {
    res += (trips[trip_num].acc_distance[i2] - trips[trip_num].acc_distance[i1+2]) * diff;
    res -= trips[trip_num].acc_penalty[i1+2] - trips[trip_num].acc_penalty[i1];
    res -= trips[trip_num].acc_penalty[i2+2] - trips[trip_num].acc_penalty[i2];
    if (res > max_margin) return false;
    res += dist(gift_prev, trips[trip_num].seq[i2+1]) * trips[trip_num].acc_weight[i1];
    if (res > max_margin) return false;
    res += dist(trips[trip_num].seq[i2+1], trips[trip_num].seq[i1+2]) * (trips[trip_num].acc_weight[i1+1] + diff);
    if (res > max_margin) return false;
    res += dist(trips[trip_num].seq[i2], trips[trip_num].seq[i1+1]) * (trips[trip_num].acc_weight[i2] + diff);
    if (res > max_margin) return false;
    res += dist(trips[trip_num].seq[i1+1], gift_next) * trips[trip_num].acc_weight[i2+1];
  }
  *real_margin = res;
  return res <= max_margin;
}

bool acceptable(double diff, double threshold, CRandomMother* _random, double prob_multiplier) {
  if (diff <= 0) return true;
  if (diff >= -threshold * 17) return false;
  return exp(diff / threshold) * prob_multiplier > _random->Random();
}

double acceptable_max_margin(double threshold, CRandomMother* _random, double prob_multiplier) {
  return threshold * log(_random->Random() / prob_multiplier + 0.000000001);
}

void save_partial_result(char* out_file, const vector<vector<int> >& seq) {
  FILE *out = fopen(out_file, "w");
  fprintf(out, "GiftId,TripId\n");
  for (int i=0; i<gift_cnt; i++) {
    for (int j=2; j<seq[i].size(); j++)
      fprintf(out, "%d,%d\n", seq[i][j-1]+1, i);
  }
  fclose(out);
}

void sa_work(const char* progress_filename, const long long iter, const double start_temp, const int thread_num, const int thread_cnt) {
  CRandomMother* _random = new CRandomMother(time(0) + thread_num * 2);

  FILE *o3;
  if (thread_num == 0)
    o3 = fopen(progress_filename, "w");

  int gift_start = (thread_num * gift_cnt) / thread_cnt;
  int gift_step = ((thread_num + 1) * gift_cnt) / thread_cnt - gift_start;
  double temp_now = start_temp;

  double start_time = getTime();
  double prv_time = start_time;

  for (long long tr=0; ; tr++) {
    if ((!(tr&((1<<20)-1)) || tr == iter)) {
      thread_iter[thread_num] = tr;
      if (thread_num == 0) new_change = -start_temp * (1 - 1.0 * tr / iter);
      temp_now = new_change;
      if (thread_num == 0 && (!(tr&((1<<27)-1)) || tr == iter)) {
        // some very inefficient operations which're called very infrequently.
        vector<unique_lock<mutex> > locks;
        for (int i=0; i<gift_cnt; i++)
          locks.emplace_back(trips[i].trip_mutex, defer_lock);
        for (int i=0; i<gift_cnt; i++)
          locks[i].lock();
        vector<vector<int> > seq(gift_cnt);
        vector<double> wei(gift_cnt);
        for (int i=0; i<gift_cnt; i++) {
          seq[i] = trips[i].seq;
          wei[i] = trips[i].weight_sum;
        }
        double gain_sum = 0.0;
        for (int i=0; i<thread_cnt; i++) {
          gain_sum += gain[i] * 6371;
        }
        for (int i=0; i<gift_cnt; i++)
          locks[i].unlock();

        double res = 0.0;
        vector<int> sanity(gift_cnt);
        for (int i=0; i<gift_cnt; i++) {
          res += calc_penalty(seq[i], wei[i]);
          for (int j=1; j+1<seq[i].size(); j++)
            sanity[seq[i][j]] ++;
        }
        for (int i=0; i<gift_cnt; i++)
          if (sanity[i] != 1) {
            printf("insane %d!\n", i);
            new_change = 0;
            break;
          }

        if (!(tr&((1<<30)-1))) {
          char tt[40]; sprintf(tt, "%s_%07lld.csv", progress_filename, thread_iter[thread_num] / 1024 / 1024);
          save_partial_result(tt, seq);
        }

        double now_time = getTime();

        printf("%.3f :", res * 6371);
        fprintf(o3, "%.3f :", res * 6371);
        for (int i=0; i<thread_cnt; i++) {
          printf(" %7lld", thread_iter[i] / 1024 / 1024);
          fprintf(o3, " %7lld", thread_iter[i] / 1024 / 1024);
        }
        printf(" [%.3f]", res * 6371 - gain_sum);
        fprintf(o3," [%.3f]", res * 6371 - gain_sum);

        printf(" %9.2fs %9.2fs\n", now_time - prv_time, now_time - start_time);
        fprintf(o3," %9.2fs %9.2fs\n", now_time - prv_time, now_time - start_time);

        prv_time = now_time;
      }
      if (temp_now == 0) break;
    }

start:
    int aaa = -1, bbb = -1, a, b;
    uint32_t rd2 = _random->BRandom();
    aaa = (rd2 >> (data_used_bit + 1)) % gift_step + gift_start;
    double prob_multiplier = 1.0;

    {
      int ord;
      if (rd2 & (1 << data_used_bit))
        ord = rd2 & ((1 << (data_used_bit - 1)) - 1);
      else
        ord = rd2 & ((1 << data_used_bit) - 1);
      bbb = db_by_dist[aaa][ord];
      if (ord >= (1 << (data_used_bit - 1)))
        prob_multiplier = 3.0;
    }

    a = cur_trip[aaa] >> 11;
    b = cur_trip[bbb] >> 11;
    if (a == b) {
      std::unique_lock<std::mutex> lock_a(trips[a].trip_mutex, std::defer_lock);
      if (!lock_a.try_lock()) goto start;
      int i1 = cur_trip[aaa]&((1<<11)-1);
      int i2 = cur_trip[bbb]&((1<<11)-1);
      int a_sz = trips[a].seq.size();
      if (a_sz <= i1+1 || trips[a].seq[i1+1] != aaa) continue;
      if (a_sz <= i2+1 || trips[a].seq[i2+1] != bbb) continue;

      if (!(rd2 & (3u<<30))) {
        int scope = 10;
        if (!(rd2 & (3u<<28))) scope = 200;
        int st = max(i1 - scope, 0);
        int fi = min(i1 + scope, a_sz - 3);
        int i2 = _random->BRandom() % (fi-st) + st;
        if (i2 >= i1) i2 ++;
      }
      int i3 = i2 + ((rd2 & (1u<<26)) >> 26);
      if ((i3 >= i1-1 && i3 <= i1+2) || (rd2 & (1u<<27))) {
        // switch
        double margin = acceptable_max_margin(temp_now, _random, prob_multiplier);
        double real_margin;
        if (check_delta_swapped_gift(a, i1, i2, margin, &real_margin)) {
          gain[thread_num] += real_margin;
          swap(trips[a].seq[i1+1], trips[a].seq[i2+1]);
          trips[a].step_distance[i1] = dist(trips[a].seq[i1], trips[a].seq[i1+1]);
          trips[a].step_distance[i1+1] = dist(trips[a].seq[i1+2], trips[a].seq[i1+1]);
          trips[a].step_distance[i2] = dist(trips[b].seq[i2], trips[a].seq[i2+1]);
          trips[a].step_distance[i2+1] = dist(trips[b].seq[i2+2], trips[a].seq[i2+1]);

          update_trip_without_distance(a, min(i1, i2), min(i1, i2));
        }
      } else {
        // insert
        i2 = i3;
        double diff = delta_moved_gift_same_trip(a, i1, i2);
        if (acceptable(diff, temp_now, _random, prob_multiplier)) {
//          printf("%d %d %f %f %f", i1, i2, diff, trips[a].cur_penalty, diff + trips[a].cur_penalty);
          gain[thread_num] += diff;
          int present = trips[a].seq[i1+1];
          if (i2 > i1) {
            for (int i=i1; i<i2-1; i++) { trips[a].seq[i+1] = trips[a].seq[i+2]; trips[a].step_distance[i+1] = trips[a].step_distance[i+2]; }
            trips[a].seq[i2] = present;
            trips[a].step_distance[i1] = dist(trips[a].seq[i1], trips[a].seq[i1+1]);
            trips[a].step_distance[i2-1] = dist(trips[a].seq[i2-1], trips[a].seq[i2]);
            trips[a].step_distance[i2] = dist(trips[a].seq[i2], trips[a].seq[i2+1]);
          } else {
            for (int i=i1; i>i2; i--) { trips[a].seq[i+1] = trips[a].seq[i]; trips[a].step_distance[i+1] = trips[a].step_distance[i]; }
            trips[a].seq[i2+1] = present;
            trips[a].step_distance[i1+1] = dist(trips[a].seq[i1+2], trips[a].seq[i1+1]);
            trips[a].step_distance[i2] = dist(trips[a].seq[i2], trips[a].seq[i2+1]);
            trips[a].step_distance[i2+1] = dist(trips[a].seq[i2+1], trips[a].seq[i2+2]);
          }
//          printf(" %f\n", trips[a].cur_penalty);
          update_trip_without_distance(a, min(i1, i2), min(i1, i2));
        }
      }
      continue;
    }

    if (!(rd2 & (127u<<25))) {
      a = _random->BRandom() % gift_cnt;
      aaa = bbb = -1;
      if (a == b) continue;
    }

    std::unique_lock<std::mutex> lock_a(trips[a].trip_mutex, std::defer_lock);
    std::unique_lock<std::mutex> lock_b(trips[b].trip_mutex, std::defer_lock);
    int x = try_lock(lock_a, lock_b);
    if (x != -1) goto start;

    if (trips[a].weight_sum < trips[b].weight_sum && (rd2&(1<<29))) {
      swap(a, b);
      swap(aaa, bbb);
    }

    int a_sz = trips[a].seq.size();
    int b_sz = trips[b].seq.size();

    if (b_sz <= 2) continue;

    int t1 = -1, t2 = -1, t3 = -1, t4 = -1;
    if (aaa > -1 && bbb > -1) {
      t1 = cur_trip[aaa]&((1<<11)-1);
      t2 = cur_trip[bbb]&((1<<11)-1);
      if (a_sz <= t1+1 || trips[a].seq[t1+1] != aaa) continue;
      if (b_sz <= t2+1 || trips[b].seq[t2+1] != bbb) continue;
    }

    vector<int> aa, bb;
    double w1, w2;
    uint32_t seed = _random->BRandom();
    int mode = (seed >> 11) % 4;
    uint32_t mm = seed >> 22;

    double delta_a = -1, delta_b = -1;
    if (mode == 0) {
      if ((mm & 3) && aaa > -1 && bbb > -1) {
        t1 += (mm >> 2) & 1;
      } else {
        uint32_t rrd = _random->BRandom();
        t2 = rrd % (b_sz - 2);
        t1 = (rrd >> 16) % (a_sz - 1);
      }
      double ww = gifts[trips[b].seq[t2+1]].weight;
      w1 = trips[a].weight_sum + ww;
      w2 = trips[b].weight_sum - ww;
      if (w1 < 1000) {
        delta_a = delta_inserted_gift(a, t1, trips[b].seq[t2+1], ww);
        delta_b = delta_removed_gift(b, t2, ww);
      } else {
        // second chance
        uint32_t rrd = _random->BRandom();
        t3 = rrd % (a_sz - 2);
        if (mm & 8) {
          for (int trial = 0; trial < 10; trial ++) {
            if (w1 - gifts[trips[a].seq[t3+1]].weight >= 1000) {
              if (t3 == a_sz - 3) break;
              t3 ++;
            }
          }
        } else {
          for (int trial = 0; trial < 10; trial ++) {
            if (w1 - gifts[trips[a].seq[t3+1]].weight >= 1000) {
              if (t3 == 0) break;
              t3 --;
            }
          }
        }
        int seq_t3 = trips[a].seq[t3+1];
        w1 -= gifts[seq_t3].weight;
        w2 += gifts[seq_t3].weight;
        if (w1 >= 1000 || w2 >= 1000) continue;

        for (int i = 0; i < (1 << data_used_bit); i ++)
          if (b == (cur_trip[db_by_dist[seq_t3][i]] >> 11)) {
            t4 = cur_trip[db_by_dist[seq_t3][i]]&((1<<11)-1);
            break;
          }
        if (t4 == -1) {
          t4 = (rrd >> 16) % (b_sz - 2);
        }

        delta_a = delta_interchanged_gift(a, t1, trips[b].seq[t2+1], t3);
        delta_b = delta_interchanged_gift(b, t4, trips[a].seq[t3+1], t2);
      }
    } else if (mode == 1) {
      if (a_sz <= 2) continue;
      if (((mm >> 2) & 3) && aaa > -1 && bbb > -1) {
      } else if (aaa > -1 && bbb > -1) {
        t2 = -1;
        if (mm & 1) {
          for (int i = 0; i < (1 << data_used_bit); i ++)
            if (b == (cur_trip[db_by_dist[aaa][i]] >> 11)) {
              t2 = cur_trip[db_by_dist[aaa][i]]&((1<<11)-1);
              break;
            }
        }
        if (t2 == -1) {
          uint32_t rrd = _random->BRandom();
          t2 = (rrd >> 16) % (b_sz - 2);
        }
      } else {
        uint32_t rrd = _random->BRandom();
        t1 = rrd % (a_sz - 2);
        t2 = (rrd >> 16) % (b_sz - 2);
      }

      w1 = trips[a].weight_sum - (gifts[trips[a].seq[t1+1]].weight - gifts[trips[b].seq[t2+1]].weight);
      w2 = trips[b].weight_sum + (gifts[trips[a].seq[t1+1]].weight - gifts[trips[b].seq[t2+1]].weight);
      if (w1 >= 1000 || w2 >= 1000) continue;

      delta_a = delta_changed_gift(a, t1, trips[b].seq[t2+1]);
      delta_b = delta_changed_gift(b, t2, trips[a].seq[t1+1]);
    } else if (mode == 2) {
      if (a_sz <= 2) continue;
      if ((mm & 12) && aaa > -1 && bbb > -1) {
        t1 += mm & 1;
        t2 += (mm >> 1) & 1;
      } else {
        uint32_t rrd = _random->BRandom();
        t1 = rrd % (a_sz - 2) + (mm & 1);
        t2 = (rrd >> 16) % (b_sz - 2) + ((mm >> 1) & 1);
      }

      if (t1 == 0 && t2 == 0) continue;

      // mix
      w1 = trips[a].weight_sum - (trips[a].acc_weight[t1] - trips[b].acc_weight[t2]);
      w2 = trips[b].weight_sum + (trips[a].acc_weight[t1] - trips[b].acc_weight[t2]);
/*      for (int trial = 0; trial < 5; trial ++) {
        if (w2 >= 1000) {
          if (t1 == a_sz) break;
          t1 ++;
          w1 -= trips[a].acc_weight[t1] - trips[a].acc_weight[t1-1];
          w2 += trips[a].acc_weight[t1] - trips[a].acc_weight[t1-1];
        }
        if (w1 >= 1000) {
          if (t2 == b_sz) break;
          t2 ++;
          w2 -= trips[b].acc_weight[t2] - trips[b].acc_weight[t2-1];
          w1 += trips[b].acc_weight[t2] - trips[b].acc_weight[t2-1];
        }
      }*/
      if (t1 == a_sz-2 && t2 == b_sz-2) continue;
      if (w1 >= 1000 || w2 >= 1000) continue;

      delta_a = delta_changed_trip(a, t1, a_sz - 3, b, t2, b_sz - 3);
      delta_b = delta_changed_trip(b, t2, b_sz - 3, a, t1, a_sz - 3);
    } else {
      if (a_sz <= 2) continue;
      if (aaa == -1 || bbb == -1) continue;

      uint32_t rrd = _random->BRandom();
      t3 = rrd % (a_sz - 2);
      int seq_t3 = trips[a].seq[t3+1];
      t4 = -1;
      for (int i = 0; i < (1 << data_used_bit); i ++)
        if (b == (cur_trip[db_by_dist[seq_t3][i]] >> 11)) {
          t4 = cur_trip[db_by_dist[seq_t3][i]]&((1<<11)-1);
          break;
        }
      if (t4 == -1) {
        t4 = (rrd >> 16) % (b_sz - 2);
      }
      if (t1 > t3) swap(t1, t3);
      if (t2 > t4) swap(t2, t4);

      w1 = trips[a].weight_sum + (trips[a].acc_weight[t3+1]-trips[a].acc_weight[t1]) - (trips[b].acc_weight[t4+1]-trips[b].acc_weight[t2]);
      w2 = trips[a].weight_sum + trips[b].weight_sum - w1;

      if (mm & 1) {
        for (int trial = 0; trial < 10; trial ++) {
          if (w1 >= 1000) {
            if (t3 == a_sz - 3) break;
            t3 ++;
            w1 += trips[a].acc_weight[t3+1] - trips[a].acc_weight[t3];
            w2 -= trips[a].acc_weight[t3+1] - trips[a].acc_weight[t3];
          }
          if (w2 >= 1000) {
            if (t4 == b_sz - 3) break;
            t4 ++;
            w2 += trips[b].acc_weight[t4+1] - trips[b].acc_weight[t4];
            w1 -= trips[b].acc_weight[t4+1] - trips[b].acc_weight[t4];
          } else if (w1 < 1000) break;
        }
      } else {
        for (int trial = 0; trial < 10; trial ++) {
          if (w1 >= 1000) {
            if (t1 == 0) break;
            t1 --;
            w1 += trips[a].acc_weight[t1+1] - trips[a].acc_weight[t1];
            w2 -= trips[a].acc_weight[t1+1] - trips[a].acc_weight[t1];
          }
          if (w2 >= 1000) {
            if (t2 == 0) break;
            t2 --;
            w2 += trips[b].acc_weight[t2+1] - trips[b].acc_weight[t2];
            w1 -= trips[b].acc_weight[t2+1] - trips[b].acc_weight[t2];
          } else if (w1 < 1000) break;
        }
      }

      if (w1 >= 1000 || w2 >= 1000) continue;

      prob_multiplier = 3.0;

      delta_a = delta_changed_trip(a, t1, t3, b, t2, t4);
      delta_b = delta_changed_trip(b, t2, t4, a, t1, t3);
    }

    double diff = delta_a + delta_b;

    if (acceptable(diff, temp_now, _random, prob_multiplier)) {
      gain[thread_num] += diff;
      trips[a].weight_sum = w1;
      trips[b].weight_sum = w2;
      if (mode == 0 && t3 == -1) {
//        printf("inserted present %d of trip %d before present %d of trip %d, gaining %.3f points.\n",
//               trips[b].seq[t2+1], b, trips[a].seq[t1+1], a, -diff * 6371);
        trips[a].seq.resize(a_sz + 1);
        trips[a].step_distance.resize(a_sz);
        for (int i=a_sz-2; i>t1; i--) {
          trips[a].seq[i+1] = trips[a].seq[i];
          trips[a].step_distance[i+1] = trips[a].step_distance[i];
        }
        trips[a].seq[a_sz] = gift_cnt;
        trips[a].seq[t1+1] = trips[b].seq[t2+1];
        trips[a].step_distance[t1] = dist(trips[a].seq[t1], trips[a].seq[t1+1]);
        trips[a].step_distance[t1+1] = dist(trips[a].seq[t1+2], trips[a].seq[t1+1]);

        for (int i=t2; i<b_sz-2; i++) {
          trips[b].seq[i+1] = trips[b].seq[i+2];
          trips[b].step_distance[i] = trips[b].step_distance[i+1];
        }
        int g1, g2;
        g1 = trips[b].seq[t2];
        g2 = trips[b].seq[t2+1];
        trips[b].step_distance[t2] = dist(g1, g2);
        trips[b].seq.resize(b_sz - 1);
        trips[b].step_distance.resize(b_sz - 2);
        update_trip_without_distance(a, 0, t1);
        update_trip_without_distance(b, 0, t2);
      } else if (mode == 0) {
//        printf("inserted present %d of trip %d before present %d of trip %d, also inserted present %d of trip %d before present %d of trip %d, gaining %.3f points.\n",
//               trips[b].seq[t2+1], b, trips[a].seq[t1+1], a, trips[a].seq[t3+1], a, trips[b].seq[t4+1], b, -diff * 6371);
        vector<int> aa, bb;
        for (int i=0; i<a_sz; i++) {
          if (i == t1+1) aa.push_back(trips[b].seq[t2+1]);
          if (i != t3+1) aa.push_back(trips[a].seq[i]);
        }
        for (int i=0; i<b_sz; i++) {
          if (i == t4+1) bb.push_back(trips[a].seq[t3+1]);
          if (i != t2+1) bb.push_back(trips[b].seq[i]);
        }
        trips[a].seq = aa;
        trips[b].seq = bb;
        update_trip(a, 0);
        update_trip(b, 0);
      } else if (mode == 1) {
//        printf("swapped present %d of trip %d and present %d of trip %d, gaining %.3f points.\n",
//               trips[a].seq[t1+1], a, trips[b].seq[t2+1], b, -diff * 6371);
        swap(trips[a].seq[t1+1], trips[b].seq[t2+1]);
        trips[a].step_distance[t1] = dist(trips[a].seq[t1], trips[a].seq[t1+1]);
        trips[a].step_distance[t1+1] = dist(trips[a].seq[t1+2], trips[a].seq[t1+1]);
        trips[b].step_distance[t2] = dist(trips[b].seq[t2], trips[b].seq[t2+1]);
        trips[b].step_distance[t2+1] = dist(trips[b].seq[t2+2], trips[b].seq[t2+1]);
        update_trip_without_distance(a, 0, t1);
        update_trip_without_distance(b, 0, t2);
      } else if (mode == 2) {
//        printf("changed trip %d from present %d to end, and trip %d from present %d to end, gaining %.3f points.\n",
//               a, trips[a].seq[t1+1], b, trips[b].seq[t2+1], -diff * 6371);
        aa.resize(t1+1);
        for (int i=0; i<t1+1; i++) { aa[i] = trips[a].seq[i]; }
        for (int i=t2+1; i<b_sz; i++) { aa.push_back(trips[b].seq[i]); }
        bb.resize(t2+1);
        for (int i=0; i<t2+1; i++) { bb[i] = trips[b].seq[i]; }
        for (int i=t1+1; i<a_sz; i++) { bb.push_back(trips[a].seq[i]); }

        int g1, g2;
        g1 = trips[a].seq[t1];
        g2 = trips[b].seq[t2+1];
        trips[a].step_distance[t1] = dist(g1, g2);
        g1 = trips[b].seq[t2];
        g2 = trips[a].seq[t1+1];
        trips[b].step_distance[t2] = dist(g1, g2);
        g1 = t1 + 1; g2 = t2 + 1;
        while (g1 <= a_sz-2 && g2 <= b_sz-2) {
          swap(trips[a].step_distance[g1], trips[b].step_distance[g2]);
          g1 ++;
          g2 ++;
        }
        while (g1 <= a_sz-2)
          trips[b].step_distance.push_back(trips[a].step_distance[g1++]);
        while (g2 <= b_sz-2)
          trips[a].step_distance.push_back(trips[b].step_distance[g2++]);

        trips[a].seq = aa;
        trips[b].seq = bb;
        trips[a].step_distance.resize(aa.size() - 1);
        trips[b].step_distance.resize(bb.size() - 1);
        update_trip_without_distance(a, 0, t1);
        update_trip_without_distance(b, 0, t2);
      } else {
//        printf("changed trip %d from present %d to present %d, and trip %d from present %d to present %d, gaining %.3f points.\n",
///               a, trips[a].seq[t1+1], trips[a].seq[t3+1], b, trips[b].seq[t2+1], trips[b].seq[t4+1], -diff * 6371);
        for (int i=0; i<t1+1; i++) { aa.push_back(trips[a].seq[i]); }
        for (int i=t2+1; i<=t4+1; i++) { aa.push_back(trips[b].seq[i]); }
        for (int i=t3+2; i<a_sz; i++) { aa.push_back(trips[a].seq[i]); }
        for (int i=0; i<t2+1; i++) { bb.push_back(trips[b].seq[i]); }
        for (int i=t1+1; i<=t3+1; i++) { bb.push_back(trips[a].seq[i]); }
        for (int i=t4+2; i<b_sz; i++) { bb.push_back(trips[b].seq[i]); }
        trips[a].seq = aa;
        trips[b].seq = bb;
        update_trip(a, 0);
        update_trip(b, 0);
      }
    }
  }

  if (thread_num == 0) fclose(o3);
}

void save_result(char* out_file, char* desc_file) {
  FILE *out = fopen(out_file, "w");
  fprintf(out, "GiftId,TripId\n");
  for (int i=0; i<gift_cnt; i++) {
    for (int j=2; j<trips[i].seq.size(); j++)
      fprintf(out, "%d,%d\n", trips[i].seq[j-1]+1, i);
  }
  fclose(out);

  out = fopen(desc_file, "w");
  fprintf(out, "GiftId,TripId,Latitude,Longitude,Weight\n");
  for (int i=0; i<gift_cnt; i++) {
    for (int j=2; j<trips[i].seq.size(); j++)
      fprintf(out, "%d,%d,%f,%f,%f\n", trips[i].seq[j-1]+1, i, gifts[trips[i].seq[j-1]].lat/PI*180, gifts[trips[i].seq[j-1]].lng/PI*180, gifts[trips[i].seq[j-1]].weight);
  }
}

void read_data() {
  FILE *in = fopen("gifts.csv", "r");
  for (int i=0; i<gift_cnt; i++) {
    int id;
    fscanf(in, "%d,%lf,%lf,%lf", &id, &gifts[i].lat, &gifts[i].lng, &gifts[i].weight);
    gifts[i].lat *= 2 * PI / 360;
    gifts[i].lng *= 2 * PI / 360;
  }
  gifts[gift_cnt].lat = PI / 2;
  gifts[gift_cnt].weight = 10;
  fclose(in);

  for (int i=0; i<=gift_cnt; i++) {
    gifts[i].x = 0.5 * cos(gifts[i].lat) * sin(gifts[i].lng);
    gifts[i].y = 0.5 * cos(gifts[i].lat) * cos(gifts[i].lng);
    gifts[i].z = 0.5 * sin(gifts[i].lat);
  }

  FILE *in2 = fopen("prv.csv", "r");
  char tmp[50]; fscanf(in2, "%s", tmp);
  for (int i=0; i<gift_cnt; i++) trips[i].seq.push_back(gift_cnt);
  for (int i=0; i<gift_cnt; i++) {
    int a, b;
    fscanf(in2, "%d,%d", &a, &b);
    trips[b].seq.push_back(a-1);
    trips[b].weight_sum += gifts[a-1].weight;
  }
  for (int i=0; i<gift_cnt; i++) trips[i].seq.push_back(gift_cnt);

  FILE *in3 = fopen("sort_by_dist", "r");
  for (int i=0; i<gift_cnt; i++) {
    char str[data_cnt*3+1];
    fscanf(in3, "%s", str);
    for (int j=0; j<(1<<data_used_bit); j++) {
      int k = (str[j*3]-'0')<<12;
      k += (str[j*3+1]-'0')<<6;
      k += str[j*3+2]-'0';
      db_by_dist[i][j] = k;
    }
  }
  fclose(in3);
}

int main(int argc, char** argv)
{
  if (argc < 7) return 0;
  printf("output = %s\n", argv[1]);
  printf("desc = %s\n", argv[2]);
  printf("progress = %s\n", argv[3]);
  printf("iter = %s\n", argv[4]);
  printf("temp = %s\n", argv[5]);
  printf("thread cnt = %s\n", argv[6]);

  long long iter;
  sscanf(argv[4], "%lld", &iter);
  double start_temp;
  sscanf(argv[5], "%lf", &start_temp);
  new_change = -start_temp;
  int thread_cnt;
  sscanf(argv[6], "%d", &thread_cnt);

  read_data();
  fprintf(stderr, "Read complete\n");

  for (int i=0; i<gift_cnt; i++) {
    update_trip(i, 0);
  }

  thread *tt = new thread[thread_cnt];
  for (int i=0; i<thread_cnt; i++)
    tt[i] = thread(sa_work, argv[3], iter, start_temp, i, thread_cnt);

  for (int i=thread_cnt-1; i>=0; i--)
    tt[i].join();

  save_result(argv[1], argv[2]);
}
