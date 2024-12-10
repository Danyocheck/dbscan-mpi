#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <mpi.h>
#include "uthash/src/uthash.h"
// #include <uthash.h>

#define EPSILON 0.05    // Пороговое расстояние
#define MIN_POINTS 5   // Минимальное количество точек для формирования кластера
#define DIMENSIONS 2

// Структура для хранения информации о пересекающихся кластерах
typedef struct {
    int cluster_id;
    int root_id;
    UT_hash_handle hh;
} ClusterMap;

typedef struct {
    double coords[DIMENSIONS]; // Координаты точки
    int cluster_id; // Идентификатор кластера
    int is_boundary; // Флаг граничной точки
    int region_id;
} Point;

// Структура для описания многомерного региона
typedef struct {
    double min_coords[DIMENSIONS]; // Минимальные координаты региона
    double max_coords[DIMENSIONS]; // Максимальные координаты региона
} Region;

// Функция для вычисления евклидова расстояния между двумя многомерными точками
double distance(Point a, Point b) {
    double dist = 0.0;
    for (int i = 0; i < DIMENSIONS; i++) {
        dist += (a.coords[i] - b.coords[i]) * (a.coords[i] - b.coords[i]);
    }
    return sqrt(dist);
}

// Функция для чтения точек из CSV файла
void read_points_from_file(const char *filename, Point **points, Region **initial_region, int *num_points) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Пропускаем заголовок
    char line[256];
    fgets(line, sizeof(line), file);

    // Выделяем начальную память для точек
    int capacity = 100;
    *points = (Point *)malloc(capacity * sizeof(Point));
    if (!*points) {
        perror("Error allocating memory for points");
        exit(EXIT_FAILURE);
    }

    *num_points = 0;
    while (fgets(line, sizeof(line), file)) {
        if (*num_points >= capacity) {
            capacity *= 2;
            *points = (Point *)realloc(*points, capacity * sizeof(Point));
            if (!*points) {
                perror("Error reallocating memory for points");
                exit(EXIT_FAILURE);
            }
        }

        sscanf(line, "%lf,%lf", &(*points)[*num_points].coords[0], &(*points)[*num_points].coords[1]);
        (*points)[*num_points].cluster_id = 0;
        (*points)[*num_points].is_boundary = 0;
        (*points)[*num_points].region_id = 0;
        (*num_points)++;
    }

    fclose(file);

    // Выделяем память для региона
    *initial_region = (Region *)malloc(sizeof(Region));
    if (!*initial_region) {
        perror("Error allocating memory for region");
        exit(EXIT_FAILURE);
    }

    // Инициализация минимальных и максимальных координат
    for (int d = 0; d < DIMENSIONS; d++) {
        (*initial_region)->min_coords[d] = DBL_MAX;
        (*initial_region)->max_coords[d] = DBL_MIN;
    }

    // Обновление минимальных и максимальных координат
    for (int i = 0; i < *num_points; i++) {
        for (int d = 0; d < DIMENSIONS; d++) {
            if ((*points)[i].coords[d] < (*initial_region)->min_coords[d]) {
                (*initial_region)->min_coords[d] = (*points)[i].coords[d];
            }
            if ((*points)[i].coords[d] > (*initial_region)->max_coords[d]) {
                (*initial_region)->max_coords[d] = (*points)[i].coords[d];
            }
        }
    }
}

// Функция для выделения региона с учетом граничной области
void get_region_data(Point **points, int *num_points, Region *initial_region, int region_id, int num_regions) {
    // Рассчитываем количество регионов по каждой оси
    int regions_per_side = (int)(sqrt(num_regions));

    Region region;

    // Инициализируем массив делений для каждой оси
    int n = (int)log2(num_regions);
    int *divisions = (int *)malloc(DIMENSIONS * sizeof(int));
    int baseDivisions = n / DIMENSIONS;  // Базовое количество делений
    int remainder = n % DIMENSIONS;     // Остаток от деления

    for (int d = 0; d < DIMENSIONS; ++d) {
        divisions[d] = 1 << baseDivisions; // Начальные деления
        if (remainder > 0) {               // Распределяем остаток
            divisions[d] *= 2;
            --remainder;
        }
    }

    // Вычисляем координаты для каждого измерения
    int region_id_tmp = region_id;
    for (int d = 0; d < DIMENSIONS; ++d) {
        int region_index = region_id_tmp % divisions[d];
        region_id_tmp /= divisions[d];
        double range = initial_region->max_coords[d] - initial_region->min_coords[d]; // Длина гиперкуба в данном измерении
        double step = range / divisions[d];        // Размер подгиперкуба
        region.min_coords[d] = initial_region->min_coords[d] + region_index * step;
        region.max_coords[d] = region.min_coords[d] + step;
    }

    free(divisions);

    // Новый массив точек, принадлежащих подрегиону
    Point *region_points = (Point *)malloc(*num_points * sizeof(Point));
    int region_point_count = 0;

    // Поиск точек в пределах региона
    for (int i = 0; i < *num_points; i++) {
        int in_region = 1;
        int is_boundary = 0;

        // Проверяем, принадлежит ли точка региону
        for (int d = 0; d < DIMENSIONS; d++) {
            if ((*points)[i].coords[d] < region.min_coords[d] - EPSILON || 
                (*points)[i].coords[d] >= region.max_coords[d] + EPSILON) {
                in_region = 0;
                break;
            }
        }

        // Если точка в регионе, добавляем ее в новый массив
        if (in_region) {
            region_points[region_point_count] = (*points)[i];

            // Помечаем точку как граничную, если она на границе региона
            for (int d = 0; d < DIMENSIONS; d++) {
                if ((*points)[i].coords[d] >= region.min_coords[d] - EPSILON &&
                    (*points)[i].coords[d] < region.min_coords[d]) {
                    is_boundary = 1;
                }
                if ((*points)[i].coords[d] >= region.max_coords[d] &&
                    (*points)[i].coords[d] < region.max_coords[d] + EPSILON) {
                    is_boundary = 1;
                }
            }

            region_points[region_point_count].cluster_id = 0;
            region_points[region_point_count].is_boundary = is_boundary;
            region_points[region_point_count].region_id = region_id;
            region_point_count++;
        }
    }

    // Освобождаем старый массив точек и обновляем указатель
    free(*points);
    *points = region_points;
    *num_points = region_point_count;
}

// Расширение кластера, начиная с данной точки
void expand_cluster(Point *points, int num_points, Point *all_points, int total_points, int index, int cluster_id, int rank) {
    points[index % num_points].cluster_id = cluster_id + rank * 1000;
    int *seeds = NULL;
    int seeds_count = 0;

    // Поиск соседних точек
    for (int i = 0; i < total_points; i++) {
        if (distance(points[index % num_points], all_points[i]) <= EPSILON) {
            seeds = (int *)realloc(seeds, (seeds_count + 1) * sizeof(int));
            seeds[seeds_count++] = i;
        }
    }

    if (seeds_count < MIN_POINTS) {
        points[index % num_points].cluster_id = -1;  // Шум
    } else {
        for (int i = 0; i < seeds_count; i++) {
            int idx = seeds[i];
            if (all_points[idx].cluster_id == 0) {
                all_points[idx].cluster_id = cluster_id + rank * 1000;
                expand_cluster(points, num_points, all_points, total_points, idx, cluster_id, rank);
            }
        }
    }
    free(seeds);
}

// Функция для проверки полного совпадения координат
bool coordinates_match(double *a, double *b) {
    for (int i = 0; i < DIMENSIONS; i++) {
        if (a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

int find_root(ClusterMap *map, int cluster_id) {
    ClusterMap *entry;
    while (1) {
        HASH_FIND_INT(map, &cluster_id, entry);
        if (entry == NULL || entry->root_id == cluster_id) {
            return cluster_id;
        }
        cluster_id = entry->root_id;
    }
}

void merge_clusters(Point *points, int num_points) {
    ClusterMap *cluster_map = NULL, *entry, *tmp;

    // Фаза 1: Построение связей между кластерами
    for (int i = 0; i < num_points; i++) {
        if (points[i].is_boundary) {
            for (int j = i + 1; j < num_points; j++) {
                if (points[i].region_id == points[j].region_id) continue;
                if (coordinates_match(points[i].coords, points[j].coords)) {
                    int cluster_a = points[i].cluster_id;
                    int cluster_b = points[j].cluster_id;

                    if (cluster_a == -1 && cluster_b == -1) continue;

                    int root_a = (cluster_a != -1) ? find_root(cluster_map, cluster_a) : -1;
                    int root_b = (cluster_b != -1) ? find_root(cluster_map, cluster_b) : -1;

                    if (root_a != -1 && root_b != -1 && root_a != root_b) {
                        // Объединить два кластера
                        int min_root = root_a < root_b ? root_a : root_b;
                        int max_root = root_a > root_b ? root_a : root_b;

                        HASH_FIND_INT(cluster_map, &max_root, entry);
                        if (entry == NULL) {
                            entry = (ClusterMap *)malloc(sizeof(ClusterMap));
                            entry->cluster_id = max_root;
                            entry->root_id = min_root;
                            HASH_ADD_INT(cluster_map, cluster_id, entry);
                        } else {
                            entry->root_id = min_root;
                        }
                    } else if (root_a == -1) {
                        points[i].cluster_id = root_b;
                    } else if (root_b == -1) {
                        points[j].cluster_id = root_a;
                    }
                }
            }
        }
    }

    // Фаза 2: Обновление всех cluster_id
    for (int i = 0; i < num_points; i++) {
        if (points[i].cluster_id != -1) {
            points[i].cluster_id = find_root(cluster_map, points[i].cluster_id);
        }
    }

    // Освобождение памяти хэш-таблицы
    HASH_ITER(hh, cluster_map, entry, tmp) {
        HASH_DEL(cluster_map, entry);
        free(entry);
    }
}


int main(int argc, char **argv) {
    int rank, size;
    Point *points = NULL;
    int num_points;
    Point *all_points = NULL;
    int total_points;
    Region *initial_region = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Имя файла с входными данными
    const char *input_file = "two_circles.csv";
    const char *output_file = "clustered_points.csv";

    double start_time, end_time;

    if (rank == 0) {
        start_time = MPI_Wtime();
        read_points_from_file(input_file, &points, &initial_region, &num_points);
        end_time = MPI_Wtime();
        printf("Time to read file: %f seconds\n", end_time - start_time);
    }

    // Рассылка общего числа точек всем процессам
    MPI_Bcast(&num_points, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Рассылка начальных координат региона всем процессам
    if (rank != 0) {
        initial_region = (Region *)malloc(sizeof(Region));
    }
    MPI_Bcast(initial_region->min_coords, DIMENSIONS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(initial_region->max_coords, DIMENSIONS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Распределение точек между процессами
    int points_per_process = num_points;
    Point *local_points = (Point *)malloc(points_per_process * sizeof(Point));

    if (rank == 0) {
        for (int proc = 1; proc < size; proc++) {
            for (int i = 0; i < points_per_process; i++) {
                MPI_Send(points[i].coords, DIMENSIONS, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
            }
        }

        for (int i = 0; i < points_per_process; i++) {
            for (int d = 0; d < DIMENSIONS; d++) {
                local_points[i].coords[d] = points[i].coords[d];
            }
        }

        // Освобождение памяти на корневом процессе
        free(points);
    } else {
        for (int i = 0; i < points_per_process; i++) {
            MPI_Recv(local_points[i].coords, DIMENSIONS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    // Ожидание пока все процессы не будут готовы к кластеризации
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        printf("start clasterization\n");
        start_time = MPI_Wtime();
    }

    // Получение данных региона для каждого процесса
    get_region_data(&local_points, &points_per_process, initial_region, rank, size);

    // Кластеризация точек
    int cluster_id = 1;
    start_time = MPI_Wtime();
    for (int i = 0; i < points_per_process; i++) {
        if (local_points[i].cluster_id == 0) {
            expand_cluster(local_points, points_per_process, local_points, points_per_process, i, cluster_id, rank);
            cluster_id++;
        }
    }


    // Сбор данных всех процессов
    int *recv_counts = NULL;
    int *displs = NULL;
    if (rank == 0) {
        recv_counts = (int *)malloc(size * sizeof(int));
        displs = (int *)malloc(size * sizeof(int));
    }

    MPI_Gather(&points_per_process, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        displs[0] = 0;  // Начальное смещение для процесса 0
        for (int i = 1; i < size; i++) {
            displs[i] = recv_counts[i - 1] + displs[i - 1];  // Смещение в элементах
        }

        all_points = (Point *)malloc((recv_counts[size - 1] + displs[size - 1]) * sizeof(Point));
    }

    // Создание пользовательского типа MPI для структуры Point
    MPI_Datatype point_type;
    int block_lengths[4] = {DIMENSIONS, 1, 1, 1};
    MPI_Aint offsets[4];
    MPI_Aint base_address;
    MPI_Get_address(&local_points[0], &base_address);
    MPI_Get_address(&local_points[0].coords, &offsets[0]);
    MPI_Get_address(&local_points[0].cluster_id, &offsets[1]);
    MPI_Get_address(&local_points[0].is_boundary, &offsets[2]);
    MPI_Get_address(&local_points[0].region_id, &offsets[3]);

    for (int i = 0; i < 4; i++) {
        offsets[i] -= base_address;
    }

    MPI_Type_create_struct(4, block_lengths, offsets, 
                           (MPI_Datatype[]){MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT}, 
                           &point_type);
    MPI_Type_commit(&point_type);

    // Сбор всех точек
    MPI_Gatherv(local_points, points_per_process, point_type,
                all_points, recv_counts, displs, point_type, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Объединение пересекающихся кластеров
        merge_clusters(all_points, recv_counts[size - 1] + displs[size - 1]);
        end_time = MPI_Wtime();
        printf("Time for clustering: %f seconds\n", end_time - start_time);
        printf("end clasteriaztion\n");

        // Запись результатов
        FILE *output = fopen(output_file, "w");
        fprintf(output, "x,y,cluster_id\n");
        for (int i = 0; i < recv_counts[size - 1] + displs[size - 1]; i++) {
            if (all_points[i].is_boundary) continue;
            fprintf(output, "%f,%f,%d\n", all_points[i].coords[0], all_points[i].coords[1], all_points[i].cluster_id);
        }
        fclose(output);

        // Освобождение памяти на корневом процессе
        free(all_points);
        free(recv_counts);
        free(displs);
    }

    // Освобождение памяти на всех процессах
    free(local_points);
    free(initial_region);

    MPI_Type_free(&point_type);
    MPI_Finalize();
    return 0;
}