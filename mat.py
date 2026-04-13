import numpy as np
import matplotlib.pyplot as plt

# ファイル読み込み（double配列）
data = np.fromfile("solver5after.data", dtype=np.float64)

# n を計算
n = int(np.sqrt(len(data)))

# 2次元に変換
image = data.reshape((n, n))

# 表示
plt.imshow(image, cmap="rainbow", origin="lower")
plt.colorbar()

# 軸スケール（Igorと同じ）
plt.xlim(0, n)
plt.ylim(0, n)

plt.show()