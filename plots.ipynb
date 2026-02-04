import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==================== PARTE 1: ENTRADAS DO SISTEMA ====================
# Parâmetros
TS, TCC, PS, PD, QMAX = 0.3, 1.0, 120.0, 80.0, 300.0
TT, HT = 1.0, 0.01
tempo = np.linspace(0, TT, int(TT/HT))

# Calcula PIN e QIN
PIN, QIN = np.zeros_like(tempo), np.zeros_like(tempo)
for i, t in enumerate(tempo):
    t_ciclo = t % TCC
    if t_ciclo <= TS:
        PIN[i] = PD + ((PS - PD) / TS) * t_ciclo
        QIN[i] = QMAX * (np.sin(np.pi * t_ciclo / TS))**2
    else:
        A, B = (PS - PD) / (TS - TCC), (PD * TS - PS * TCC) / (TS - TCC)
        PIN[i] = A * t_ciclo + B
        QIN[i] = 0.0

# Plota entradas
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))

# Primeiro gráfico
ax1.plot(tempo, PIN, 'b-', linewidth=2)
ax1.set_xlabel('Tempo (s)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Pressão (mmHg)', fontsize=12, fontweight='bold', color='b')
ax1.tick_params(axis='y', labelcolor='b')
ax1.set_ylim([75, 125])
ax1.grid(True, alpha=0.3)

# Segundo gráfico
ax2.plot(tempo, QIN, color='orange', linewidth=2, linestyle='--')
ax2.set_xlabel('Tempo (s)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Vazão (mL/s)', fontsize=12, fontweight='bold', color='orange')
ax2.tick_params(axis='y', labelcolor='orange')
ax2.set_ylim([-10, 320])
ax2.grid(True, alpha=0.3)

# Título principal da figura
fig.suptitle('Entradas do Sistema (Vetores b1 e b2)', fontsize=14, fontweight='bold', y=0.98)

plt.tight_layout()
plt.show()

# ==================== PARTE 2: RESULTADOS DA SIMULAÇÃO ====================

print("\nresultados da simulação!")

# Lê dados
df = pd.read_csv('resultados_simulacao.csv')
df.columns = df.columns.str.strip()

# Plota cada trecho
for idx in [16, 30, 48, 49]:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    
    # Pressão
    ax1.plot(df['tempo'], df[f'P_{idx}'], 'b-', linewidth=1.5)
    ax1.set_xlabel('Tempo (s)', fontsize=11)
    ax1.set_ylabel('Pressão (Pa)', fontsize=11)
    ax1.set_title(f'Trecho {idx} - Pressão', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Vazão
    ax2.plot(df['tempo'], df[f'Q_{idx}'], 'r-', linewidth=1.5)
    ax2.set_xlabel('Tempo (s)', fontsize=11)
    ax2.set_ylabel('Vazão (m³/s)', fontsize=11)
    ax2.set_title(f'Trecho {idx} - Vazão', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    plt.close()
