import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==================== ANÁLISE DOS RESULTADOS DA SIMULAÇÃO ====================

print("\nResultados da simulação!")

# Lê dados
df = pd.read_csv('resultados_simulacao.csv')
df.columns = df.columns.str.strip()

# Lista dos trechos de interesse
trechos = [16, 19, 20, 23, 24, 30, 48, 49, 52, 58]

# Extrai PIN e QIN (são as entradas do sistema)
PIN = df['PIN']
QIN = df['QIN']

# Plota PIN (Pressão de entrada) e QIN (Vazão de entrada)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))

ax1.plot(df['tempo'], PIN, 'b-', linewidth=2)
ax1.set_xlabel('Tempo (s)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Pressão (mmHg)', fontsize=12, fontweight='bold', color='b')
ax1.tick_params(axis='y', labelcolor='b')
ax1.set_ylim([75, 125])
ax1.set_xlim([0, 10])
ax1.grid(True, alpha=0.3)
ax1.set_title('Pressão de Entrada (PIN)', fontsize=12, fontweight='bold')

ax2.plot(df['tempo'], QIN, color='orange', linewidth=2, linestyle='--')
ax2.set_xlabel('Tempo (s)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Vazão (m³)', fontsize=12, fontweight='bold', color='orange')
ax2.tick_params(axis='y', labelcolor='orange')
ax2.set_ylim([-10, 320])
ax2.set_xlim([0, 10])
ax2.grid(True, alpha=0.3)
ax2.set_title('Vazão de Entrada (QIN)', fontsize=12, fontweight='bold')

plt.suptitle('Entradas do Sistema', fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()
plt.show()

# Plota cada trecho individualmente
print("\nGerando gráficos individuais...")
for idx in trechos:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    # Pressão
    ax1.plot(df['tempo'], df[f'P_{idx}'], 'b-', linewidth=1.5)
    ax1.set_xlabel('Tempo (s)', fontsize=11)
    ax1.set_ylabel('Pressão (mmHg)', fontsize=11)
    ax1.set_title(f'Trecho {idx} - Pressão', fontsize=12, fontweight='bold')
    ax1.set_xlim([0, 10])
    ax1.set_ylim([-10, 160])
    ax1.grid(True, alpha=0.3)

    # Vazão
    ax2.plot(df['tempo'], df[f'Q_{idx}'], 'r-', linewidth=1.5)
    ax2.set_xlabel('Tempo (s)', fontsize=11)
    ax2.set_ylabel('Vazão (M³/s)', fontsize=11)
    ax2.set_title(f'Trecho {idx} - Vazão', fontsize=12, fontweight='bold')
    ax2.set_xlim([0, 10])
    ax2.set_ylim([-10, 320])
    ax2.grid(True, alpha=0.3)

    plt.suptitle(f'Resultados do Trecho {idx}', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.show()

print("\nAnálise concluída!")
