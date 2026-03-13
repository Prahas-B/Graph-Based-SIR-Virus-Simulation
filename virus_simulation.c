#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#define MAX_NODES 100
#define QSIZE 500

#define INPUT_FILE "graph_input.txt"
#define OUTPUT_FILE "graph_output.txt"


#define SUSCEPTIBLE 0
#define INFECTED 1
#define RECOVERED 2

typedef struct gnode
{
    int info;
    struct gnode *next;
} GNODE;

typedef struct
{
    int id;
    int state;
    float infection_severity;  
    int degree;
} PERSON;

typedef struct
{
    int id;
    int degree;
} HeapNode;

GNODE *adlist[MAX_NODES];
PERSON population[MAX_NODES];
int n;
float transmission_prob;
float decay = 0.1f;
float threshold = 0.15f;

int queue[QSIZE];
int front = -1, rear = -1;

int susceptible_count, infected_count, recovered_count;

HeapNode heap[MAX_NODES];
int heap_count;

void load_graph_from_file() {
    FILE *fp = fopen(INPUT_FILE, "r");
    if (!fp) {
        printf("\nError: Could not open %s\n", INPUT_FILE);
        return;
    }

    printf("\nLoading contact network from %s...\n", INPUT_FILE);

    if (fscanf(fp, "%d", &n) != 1) {
        printf("Error: Could not read number of nodes.\n");
        fclose(fp);
        return;
    }

    if (n < 1 || n > MAX_NODES) {
        printf("Error: n must be between 1 and %d\n", MAX_NODES);
        fclose(fp);
        return;
    }

    for (int i = 1; i <= n; i++) {
        adlist[i] = NULL;
        population[i].degree = 0;
        population[i].state = SUSCEPTIBLE;
        population[i].infection_severity = 0.0f;
        population[i].id = i;
    }

    int i, j;

    while (fscanf(fp, "%d %d", &i, &j) == 2) {
        if (i == 0 && j == 0) break;

        if (i < 1 || i > n || j < 1 || j > n) {
            printf("Invalid edge %d %d skipped\n", i, j);
            continue;
        }

        if (i == j) continue;

        int exists = 0;
        for (GNODE *check = adlist[i]; check != NULL; check = check->next) {
            if (check->info == j) exists = 1;
        }
        if (exists) continue;

        GNODE *a = malloc(sizeof(GNODE));
        a->info = j; a->next = adlist[i];
        adlist[i] = a;
        population[i].degree++;

        GNODE *b = malloc(sizeof(GNODE));
        b->info = i; b->next = adlist[j];
        adlist[j] = b;
        population[j].degree++;
    }

    fclose(fp);
    printf("Graph loaded successfully! (n = %d)\n", n);
}


void save_graph_to_file() {
    FILE *fp = fopen(OUTPUT_FILE, "w");
    if (!fp) {
        printf("\nError: Could not open %s for writing\n", OUTPUT_FILE);
        return;
    }

    fprintf(fp, "===== CONTACT NETWORK =====\n");
    for (int i = 1; i <= n; i++) {
        fprintf(fp, "Person %d [Degree %d]: ", i, population[i].degree);

        GNODE *temp = adlist[i];
        if (!temp) {
            fprintf(fp, "(No contacts)\n");
            continue;
        }
        while (temp) {
            fprintf(fp, "%d ", temp->info);
            temp = temp->next;
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\n===== POPULATION STATES =====\n");
    for (int i = 1; i <= n; i++) {
        char *state =
            (population[i].state == SUSCEPTIBLE) ? "S" :
            (population[i].state == INFECTED) ? "I" : "R";

        fprintf(fp, "Person %d: %s (Severity %.2f)\n",
                i, state, population[i].infection_severity);
    }

    fclose(fp);

    printf("\nGraph saved to %s successfully!\n", OUTPUT_FILE);
}


void creategraph();
void display();
void initialize_population();
void reset_population();
void seed_infection(int);
void run_simulation();
void simulate_day();
void recover_nodes();
void transmit_infection();
void update_timers();
void display_statistics(int);
void calculate_statistics();
int get_random_node();
void find_super_spreaders();
void heapify_degrees();
void adjust_heap(int);

void enque(int);
int deque();
int isqempty();
int isqfull();

int main()
{
    int choice, k, i;
    int graph_created = 0;

    srand((unsigned)time(NULL));

    printf("\n========== SIR MODEL CONTAGION SIMULATOR ==========\n");
    
    for(i = 1; i <= n; i++)
        adlist[i] = NULL;

    initialize_population();
    reset_population(); 

    transmission_prob = 0.3f;

    while(1)
    {
        printf("\n========== MENU ==========");
        printf("\n1. Create Contact Network (Graph)");
        printf("\n2. Display Contact Network");
        printf("\n3. Set Simulation Parameters");
        printf("\n4. Seed Infection (Patient Zero)");
        printf("\n5. Find Super-Spreaders (Top 10%%)");
        printf("\n6. Run Simulation");
        printf("\n7. Reset Simulation (Keep Graph)");
        printf("\n8. Exit");
        printf("\n==========================");
        printf("\nEnter your choice: ");
        if (scanf("%d", &choice) != 1) { printf("Invalid input.\n"); break; }

        switch(choice)
        {
            case 1:
                printf("\nLoad graph from file or enter manually?");
                printf("\n1. Load from graph_input.txt");
                printf("\n2. Enter manually");
                printf("\nChoice: ");

            int x;
            scanf("%d", &x);

            if (x == 1)
                load_graph_from_file();
            else
                creategraph();

            graph_created = 1;
                break;

            case 2:
                if(!graph_created)
                {
                    printf("\nError: Please create the contact network first!\n");
                    break;
                }
                display();
                save_graph_to_file();
                break;
            case 3:
                printf("Enter transmission probability (0.0 - 1.0): ");
                if (scanf("%f", &transmission_prob) != 1) transmission_prob = 0.3f;
                if(transmission_prob < 0.0 || transmission_prob > 1.0)
                {
                    printf("Invalid probability! Setting to default 0.3\n");
                    transmission_prob = 0.3f;
                }

                printf("\nParameters set successfully!");
                printf("\nTransmission Probability: %.2f\n", transmission_prob);
                break;
            case 4:
                if(!graph_created)
                {
                    printf("\nError: Please create the contact network first!\n");
                    break;
                }
                printf("Enter number of initial infected people (Patient Zero): ");
                if (scanf("%d", &k) != 1) { printf("Invalid input.\n"); break; }
                if(k < 1 || k > n)
                {
                    printf("Error: Number must be between 1 and %d\n", n);
                    break;
                }
                seed_infection(k);
                break;
            case 5:
                if(!graph_created)
                {
                    printf("\nError: Please create the contact network first!\n");
                    break;
                }
                find_super_spreaders();
                break;
            case 6:
                if(!graph_created)
                {
                    printf("\nError: Please create the contact network first!\n");
                    break;
                }
                run_simulation();
                break;
            case 7:
                reset_population();
                printf("\nSimulation reset! Graph structure preserved.\n");
                break;
            case 8:
                printf("\nExiting simulator. Stay safe!\n");
                return 0;
            default:
                printf("Invalid choice! Please select 1-8.\n");
        }
    }

    return 0;
}

void creategraph()
{
    int i, j;
    GNODE *temp1, *temp2, *check;
    int edge_exists;

    printf("\n--- Creating Contact Network ---\n");
    printf("Enter contacts as pairs (source destination)\n");
    printf("Example: 1 2 means person 1 has contact with person 2\n");
    printf("Enter 0 0 to stop\n\n");

    while(1)
    {
        printf("Enter source and destination: ");
        if (scanf("%d %d", &i, &j) != 2) { printf("Invalid input\n"); continue; }

        if(i == 0 && j == 0) break;

        if(i < 1 || i > n || j < 1 || j > n)
        {
            printf("Invalid node! Nodes must be between 1 and %d\n", n);
            continue;
        }

        if(i == j)
        {
            printf("Self-loops not allowed! Person cannot contact themselves.\n");
            continue;
        }

        edge_exists = 0;
        for(check = adlist[i]; check != NULL; check = check->next)
        {
            if(check->info == j)
            {
                edge_exists = 1;
                break;
            }
        }

        if(edge_exists)
        {
            printf("Edge already exists\n");
            continue;
        }

        temp1 = (GNODE*)malloc(sizeof(GNODE));
        temp1->info = j;
        temp1->next = adlist[i];
        adlist[i] = temp1;
        population[i].degree++;

        temp2 = (GNODE*)malloc(sizeof(GNODE));
        temp2->info = i;
        temp2->next = adlist[j];
        adlist[j] = temp2;
        population[j].degree++;
    }

    printf("\nContact network created successfully!\n");
}

void display()
{
    int p;
    GNODE *temp;
    char *state_str;

    printf("\n========== CONTACT NETWORK ==========\n");
    for(p = 1; p <= n; p++)
    {
        if(population[p].state == SUSCEPTIBLE)
            state_str = "SUSCEPTIBLE";
        else if(population[p].state == INFECTED)
            state_str = "INFECTED";
        else
            state_str = "RECOVERED";

        printf("Person %d [%s, Contacts: %d]: ", p, state_str, population[p].degree);

        if(adlist[p] == NULL)
        {
            printf("(No contacts)");
        }
        else
        {
            for(temp = adlist[p]; temp != NULL; temp = temp->next)
            {
                printf("%d ", temp->info);
            }
        }
        printf("\n");
    }
    printf("=====================================\n");
}

void initialize_population()
{
    int i;
    for(i = 1; i <= n; i++)
    {
        population[i].id = i;
        population[i].state = SUSCEPTIBLE;
        population[i].infection_severity = 0.0f;
        population[i].degree = 0;
    }
}

void reset_population()
{
    int i;
    for(i = 1; i <= n; i++)
    {
        population[i].state = SUSCEPTIBLE;
        population[i].infection_severity = 0.0f;
    }
    susceptible_count = n;
    infected_count = 0;
    recovered_count = 0;
}

void seed_infection(int k)
{
    int i, node, infected = 0;
    int attempts = 0;
    int max_attempts = n * 10;

    if(k > n)
    {
        printf("Error: Cannot infect more than %d people!\n", n);
        return;
    }

    int already_infected = 0;
    for(i = 1; i <= n; i++)
    {
        if(population[i].state == INFECTED)
            already_infected++;
    }

    if(already_infected + k > n)
    {
        printf("Error: Only %d susceptible people remaining!\n", n - already_infected);
        return;
    }

    printf("\n--- Seeding Infection ---\n");

    while(infected < k && attempts < max_attempts)
    {
        node = get_random_node();
        attempts++;

        if(population[node].state == SUSCEPTIBLE)
        {
            population[node].state = INFECTED;
            population[node].infection_severity = threshold + ((float)rand() / RAND_MAX) * (1.0f - threshold);
            infected++;
            printf("Person %d is now infected (Patient Zero)!\n", node);
        }
    }

    if(infected < k)
    {
        printf("Warning: Could only infect %d people after %d attempts\n", infected, attempts);
    }

    calculate_statistics();
    printf("\nCurrent Statistics:\n");
    printf("Susceptible: %d, Infected: %d, Recovered: %d\n",susceptible_count, infected_count, recovered_count);
}

void run_simulation()
{
    int day = 0;

    if(infected_count == 0)
    {
        printf("\nError: No infected people! Please seed infection first.\n");
        return;
    }

    printf("\n========== SIMULATION STARTED ==========\n");
    printf("Transmission Probability: %.2f\n", transmission_prob);
    printf("========================================\n\n");

    display_statistics(day);

    while(infected_count > 0)
    {
        day++;
        printf("\n--- DAY %d ---\n", day);
        simulate_day();
        display_statistics(day);

        if(day > 1000)
        {
            printf("\nSimulation stopped: Exceeded 1000 days (possible infinite loop)\n");
            break;
        }
    }

    printf("\n========== SIMULATION COMPLETE ==========\n");
    printf("Total days: %d\n", day);
    printf("Final - Susceptible: %d, Recovered: %d\n", susceptible_count, recovered_count);
    printf("=========================================\n");
}

void simulate_day()
{
    recover_nodes();
    transmit_infection();
    update_timers();
    calculate_statistics();
}

void recover_nodes()
{
    int i;
    for(i = 1; i <= n; i++)
    {
        if(population[i].state == INFECTED &&
           population[i].infection_severity <= threshold)
        {
            population[i].state = RECOVERED;
            printf("Person %d has RECOVERED\n", i);
        }
    }
}

void transmit_infection()
{
    int i, current, neighbor;
    GNODE *temp;
    int newly_infected[MAX_NODES];
    int new_count = 0;
    bool to_infect[MAX_NODES] = {false}; 

    front = -1;
    rear = -1;

    for(i = 1; i <= n; i++)
    {
        if(population[i].state == INFECTED)
        {
            enque(i);
        }
    }

    if(isqempty())
        return;

    while(!isqempty())
    {
        current = deque();
        if(current == -1) break;

        for(temp = adlist[current]; temp != NULL; temp = temp->next)
        {
            neighbor = temp->info;

            if(population[neighbor].state == SUSCEPTIBLE && !to_infect[neighbor])
            {
                float roll = (float)rand() / RAND_MAX;

                if(roll < transmission_prob)
                {
                    to_infect[neighbor] = true;
                    newly_infected[new_count++] = neighbor;
                }
            }
        }
    }

    for(i = 0; i < new_count; i++)
    {
        int id = newly_infected[i];
        population[id].state = INFECTED;
        population[id].infection_severity = threshold + ((float)rand() / RAND_MAX) * (1.0f - threshold);
        printf("Person %d got INFECTED\n", id);
    }
}

void update_timers()
{
    int i;
    for(i = 1; i <= n; i++)
    {
        if(population[i].state == INFECTED)
        {
            population[i].infection_severity -= decay;
            if(population[i].infection_severity < 0.0f)
                population[i].infection_severity = 0.0f;
        }
    }
}

void display_statistics(int day)
{
    printf("\n[Day %d] Susceptible: %d | Infected: %d | Recovered: %d\n", day, susceptible_count, infected_count, recovered_count);
}

void calculate_statistics()
{
    int i;
    susceptible_count = 0;
    infected_count = 0;
    recovered_count = 0;

    for(i = 1; i <= n; i++)
    {
        if(population[i].state == SUSCEPTIBLE) susceptible_count++;
        else if(population[i].state == INFECTED) infected_count++;
        else if(population[i].state == RECOVERED) recovered_count++;
    }
}

int get_random_node()
{
    return (rand() % n) + 1;
}

void find_super_spreaders()
{
    int i, k, top_count;

    printf("\n========== SUPER-SPREADER ANALYSIS ==========\n");

    int has_edges = 0;
    for(i = 1; i <= n; i++)
    {
        if(population[i].degree > 0)
        {
            has_edges = 1;
            break;
        }
    }

    if(!has_edges)
    {
        printf("No contacts in the network! All nodes are isolated.\n");
        printf("=============================================\n");
        return;
    }

    heap_count = 0;
    for(i = 1; i <= n; i++)
    {
        heap[heap_count].degree = population[i].degree;
        heap[heap_count].id = population[i].id;
        heap_count++;
    }

    heapify_degrees();

    top_count = (n * 10) / 100;
    if(top_count < 1) top_count = 1;

    printf("Top %d people with most contacts (Top 10%%):\n", top_count);
    printf("---------------------------------------------\n");

    for(k = 0; k < top_count && k < heap_count; k++)
    {
        int max_degree = heap[0].degree;
        int max_id = heap[0].id;

        if(max_degree == 0)
        {
            printf("\nRemaining people have no contacts.\n");
            break;
        }

        printf("%d. Person %d - %d contacts\n", k+1, max_id, max_degree);

        heap[0] = heap[heap_count-1];
        heap_count--;
        adjust_heap(heap_count);
    }

    printf("=============================================\n");
}

void heapify_degrees()
{
    int k, i, j;
    for(k = 1; k < heap_count; k++)
    {
        i = k;
        HeapNode key = heap[i];
        j = (i - 1) / 2;

        while(i > 0 && heap[j].degree < key.degree)
        {
            heap[i] = heap[j];
            i = j;
            j = (i - 1) / 2;
        }
        heap[i] = key;
    }
}

void adjust_heap(int count)
{
    int i, j;

    j = 0;
    i = 2 * j + 1;
    HeapNode key = heap[j];

    while(i < count)
    {
        if(i + 1 < count)
        {
            if(heap[i+1].degree > heap[i].degree)
                i++;
        }

        if(heap[i].degree > key.degree)
        {
            heap[j] = heap[i];
            j = i;
            i = 2 * j + 1;
        }
        else
            break;
    }
    heap[j] = key;
}

void enque(int x)
{
    if(isqfull())
    {
        printf("Queue Full!\n");
        return;
    }

    if(front == -1) front = 0;
    rear++;
    queue[rear] = x;
}

int deque()
{
    int x;
    if(isqempty())
    {
        printf("Queue Empty!\n");
        return -1;
    }
    x = queue[front];
    if(front == rear)
    {
        front = rear = -1;
    }
    else
    {
        front++;
    }
    return x;
}

int isqempty()
{
    if((front == -1 && rear == -1) || (front > rear))
        return 1;
    return 0;
}

int isqfull()
{
    if(rear == QSIZE - 1)
        return 1;
    return 0;
}
